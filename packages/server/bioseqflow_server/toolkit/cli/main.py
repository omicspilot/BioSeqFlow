from __future__ import annotations

"""Command-line interface for BioSeqFlow."""

from pathlib import Path

import click

from bioseqflow_server.toolkit import __version__
from bioseqflow_server.toolkit.core.config import Config
from bioseqflow_server.toolkit.preprocessing.deduplication import DuplicateRemover
from bioseqflow_server.toolkit.preprocessing.filtering import QualityFilter
from bioseqflow_server.toolkit.preprocessing.trimming import AdapterTrimmer, QualityTrimmer
from bioseqflow_server.toolkit.qc.fastqc import FastQCRunner
from bioseqflow_server.toolkit.qc.metrics import QualityMetrics
from bioseqflow_server.toolkit.qc.parsers import FastQCParser, MultiQCParser
from bioseqflow_server.toolkit.utils.io import parse_sample_sheet
from bioseqflow_server.toolkit.utils.paired_end import (
    calculate_insert_size_distribution,
    check_read_orientation,
    validate_paired_files,
)
from bioseqflow_server.toolkit.utils.parallel import process_samples_parallel


def echo_header(text: str) -> None:
    """Print a colorful header."""
    click.echo("\n" + click.style("=" * 70, fg="cyan"))
    click.echo(click.style(f"  {text}", fg="cyan", bold=True))
    click.echo(click.style("=" * 70, fg="cyan"))


def echo_success(text: str) -> None:
    """Print success message."""
    click.echo(click.style(f"[SUCCESS] {text}", fg="green", bold=True))


def echo_warning(text: str) -> None:
    """Print warning message."""
    click.echo(click.style(f"[WARNING] {text}", fg="yellow", bold=True))


def echo_error(text: str) -> None:
    """Print error message."""
    click.echo(click.style(f"[ERROR] {text}", fg="red", bold=True))


def echo_info(text: str) -> None:
    """Print info message."""
    click.echo(click.style(f"  {text}", fg="blue"))


@click.group()
@click.version_option(version=__version__)
def cli() -> None:
    """
    BioSeqFlow: Sequencing Quality Control and Preprocessing Platform

    A comprehensive toolkit for NGS data analysis with:
      - Paired-end validation and processing
      - Fuzzy adapter matching (30-50% better detection)
      - QC results parsing (FastQC/MultiQC)
      - Advanced parallel processing

    Documentation: https://omicspilot.com/projects/bioseqflow
    Examples: See examples/ directory for detailed workflows
    """
    pass


@cli.command()
@click.option(
    "--input",
    "-i",
    "input_file",
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ file",
)
@click.option(
    "--output",
    "-o",
    "output_dir",
    required=True,
    type=click.Path(),
    help="Output directory",
)
@click.option("--threads", "-t", default=1, help="Number of threads")
@click.option("--quiet", "-q", is_flag=True, help="Suppress output")
def qc(input_file: str, output_dir: str, threads: int, quiet: bool) -> None:
    """
    Run quality control analysis with FastQC.

    Examples:
      bioseqflow qc -i sample.fastq.gz -o results/ --threads 4
    """
    if not quiet:
        echo_header("Quality Control Analysis")
        echo_info(f"Input: {input_file}")
        echo_info(f"Threads: {threads}")

    runner = FastQCRunner()
    results = runner.run(input_file, output_dir, threads=threads, quiet=quiet)

    if not quiet:
        click.echo()
        if "basic_statistics" in results.get("modules", {}):
            stats = results["modules"]["basic_statistics"]
            click.echo(click.style("Basic Statistics:", fg="cyan", bold=True))
            click.echo(f"  Total sequences: {stats.get('Total Sequences', 'N/A'):,}")
            click.echo(f"  Sequence length: {stats.get('Sequence length', 'N/A')}")
            click.echo(f"  GC content: {stats.get('%GC', 'N/A')}%")

        # Calculate composite score
        metrics = QualityMetrics()
        score = metrics.calculate_composite_score(results)
        click.echo()
        click.echo(
            "Composite Quality Score: " + click.style(f"{score:.1f}/100", fg="yellow", bold=True)
        )

        if score >= 90:
            echo_success("Quality: Excellent")
        elif score >= 70:
            echo_success("Quality: Good")
        elif score >= 50:
            echo_warning("Quality: Needs Review")
        else:
            echo_error("Quality: Requires Attention")

        click.echo()
        echo_success(f"Results saved to: {output_dir}")


@cli.command()
@click.option(
    "--input",
    "-i",
    "input_file",
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ file",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    required=True,
    type=click.Path(),
    help="Output FASTQ file",
)
@click.option("--adapter", "-a", required=True, help="Adapter sequence to trim")
@click.option(
    "--fuzzy/--exact",
    default=False,
    help="Use fuzzy matching (Smith-Waterman) for 30-50% better detection",
)
@click.option(
    "--error-rate",
    "-e",
    default=0.15,
    type=float,
    help="Maximum error rate for fuzzy matching (default: 0.15)",
)
@click.option("--quiet", "-q", is_flag=True, help="Suppress output")
def trim(
    input_file: str,
    output_file: str,
    adapter: str,
    fuzzy: bool,
    error_rate: float,
    quiet: bool,
) -> None:
    """
    Trim adapter sequences from reads.

    Examples:
      # Exact matching (fast)
      bioseqflow trim -i input.fastq -o output.fastq -a AGATCGGAAGAGC

      # Fuzzy matching (30-50% better detection)
      bioseqflow trim -i input.fastq -o output.fastq -a AGATCGGAAGAGC --fuzzy
    """
    if not quiet:
        echo_header("Adapter Trimming")
        echo_info(f"Input: {input_file}")
        echo_info(f"Adapter: {adapter}")
        echo_info(f"Mode: {'Fuzzy (Smith-Waterman)' if fuzzy else 'Exact matching'}")
        if fuzzy:
            echo_info(f"Error rate: {error_rate:.0%}")

    trimmer = AdapterTrimmer()
    stats = trimmer.trim(
        input_file,
        output_file,
        adapter,
        fuzzy_match=fuzzy,
        max_error_rate=error_rate,
    )

    if not quiet:
        click.echo()
        click.echo(click.style("Results:", fg="cyan", bold=True))
        click.echo(f"  Total reads: {stats.get('total_reads', 0):,}")
        click.echo(f"  Trimmed reads: {stats.get('trimmed_reads', 0):,}")
        click.echo(
            "  Trim rate: "
            + click.style(f"{stats.get('percent_trimmed', 0):.1f}%", fg="yellow", bold=True)
        )
        click.echo(f"  Bases removed: {stats.get('total_bp_removed', 0):,} bp")

        echo_success(f"Output saved to: {output_file}")


@cli.command()
@click.option(
    "--input",
    "-i",
    "input_file",
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ file",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    required=True,
    type=click.Path(),
    help="Output FASTQ file",
)
@click.option("--min-quality", "-q", default=20, help="Minimum mean quality score")
@click.option("--min-length", "-l", default=50, help="Minimum read length")
@click.option("--max-n", default=5.0, help="Maximum N content percentage")
@click.option("--quiet", is_flag=True, help="Suppress output")
def filter(
    input_file: str,
    output_file: str,
    min_quality: int,
    min_length: int,
    max_n: float,
    quiet: bool,
) -> None:
    """Filter reads by quality and length."""
    if not quiet:
        click.echo(f"Filtering reads from {input_file}...")

    filter_module = QualityFilter()
    stats = filter_module.filter(
        input_file, output_file, min_quality=min_quality, min_length=min_length, max_n_content=max_n
    )

    if not quiet:
        click.echo("\n=== Filtering Results ===")
        click.echo(f"Total Reads: {stats.get('total_reads', 'N/A')}")
        click.echo(f"Passed Reads: {stats.get('passed_reads', 'N/A')}")
        click.echo(f"Failed Quality: {stats.get('failed_quality', 'N/A')}")
        click.echo(f"Failed Length: {stats.get('failed_length', 'N/A')}")
        click.echo(f"Failed N Content: {stats.get('failed_n_content', 'N/A')}")
        click.echo(f"Percent Passed: {stats.get('percent_passed', 0):.2f}%")
        click.echo(f"\nOutput saved to: {output_file}")


@cli.command()
@click.option(
    "--input",
    "-i",
    "input_file",
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ file",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    required=True,
    type=click.Path(),
    help="Output FASTQ file",
)
@click.option(
    "--method",
    "-m",
    type=click.Choice(["exact", "prefix"]),
    default="exact",
    help="Deduplication method",
)
@click.option("--quiet", is_flag=True, help="Suppress output")
def deduplicate(input_file: str, output_file: str, method: str, quiet: bool) -> None:
    """Remove duplicate reads."""
    if not quiet:
        click.echo(f"Removing duplicates from {input_file}...")

    dedup = DuplicateRemover()
    stats = dedup.remove_duplicates(input_file, output_file, method=method)

    if not quiet:
        click.echo("\n=== Deduplication Results ===")
        click.echo(f"Total Reads: {stats.get('total_reads', 'N/A')}")
        click.echo(f"Unique Reads: {stats.get('unique_reads', 'N/A')}")
        click.echo(f"Duplicates: {stats.get('duplicates', 'N/A')}")
        click.echo(f"Duplication Rate: {stats.get('duplication_rate', 0):.2f}%")
        click.echo(f"\nOutput saved to: {output_file}")


@cli.command()
@click.option(
    "--input",
    "-i",
    "input_file",
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ file",
)
@click.option(
    "--output",
    "-o",
    "output_dir",
    required=True,
    type=click.Path(),
    help="Output directory",
)
@click.option("--adapter", "-a", help="Adapter sequence (optional)")
@click.option("--min-quality", "-q", default=20, help="Minimum quality score")
@click.option("--min-length", "-l", default=50, help="Minimum read length")
@click.option("--threads", "-t", default=1, help="Number of threads")
@click.option("--keep-intermediate", is_flag=True, help="Keep intermediate files")
@click.option("--quiet", is_flag=True, help="Suppress output")
def pipeline(
    input_file: str,
    output_dir: str,
    adapter: str | None,
    min_quality: int,
    min_length: int,
    threads: int,
    keep_intermediate: bool,
    quiet: bool,
) -> None:
    """Run complete QC and preprocessing pipeline."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    input_path = Path(input_file)
    base_name = input_path.stem.replace(".fastq", "").replace(".fq", "")

    if not quiet:
        click.echo("=" * 60)
        click.echo("BioSeqFlow Pipeline")
        click.echo("=" * 60)

    # Step 1: Initial QC
    if not quiet:
        click.echo("\n[Step 1/5] Running initial QC...")
    qc_dir = output_path / "qc_initial"
    runner = FastQCRunner()
    initial_qc = runner.run(input_file, str(qc_dir), threads=threads, quiet=True)

    # Step 2: Adapter trimming (if adapter provided)
    current_file = input_file
    if adapter:
        if not quiet:
            click.echo("[Step 2/5] Trimming adapters...")
        trimmed_file = output_path / f"{base_name}_trimmed.fastq.gz"
        trimmer = AdapterTrimmer()
        trim_stats = trimmer.trim(current_file, str(trimmed_file), adapter)
        current_file = str(trimmed_file)
    else:
        if not quiet:
            click.echo("[Step 2/5] Skipping adapter trimming (no adapter provided)")
        trim_stats = {}

    # Step 3: Quality trimming
    if not quiet:
        click.echo("[Step 3/5] Trimming low-quality bases...")
    qtrimmed_file = output_path / f"{base_name}_qtrimmed.fastq.gz"
    qtrimmer = QualityTrimmer()
    qtrim_stats = qtrimmer.trim(current_file, str(qtrimmed_file), min_quality=min_quality)
    current_file = str(qtrimmed_file)

    # Step 4: Filtering
    if not quiet:
        click.echo("[Step 4/5] Filtering reads...")
    filtered_file = output_path / f"{base_name}_filtered.fastq.gz"
    filter_module = QualityFilter()
    filter_stats = filter_module.filter(
        current_file, str(filtered_file), min_quality=min_quality, min_length=min_length
    )
    current_file = str(filtered_file)

    # Step 5: Post-QC
    if not quiet:
        click.echo("[Step 5/5] Running post-processing QC...")
    post_qc_dir = output_path / "qc_final"
    final_qc = runner.run(current_file, str(post_qc_dir), threads=threads, quiet=True)

    # Generate summary
    if not quiet:
        click.echo("\n" + "=" * 60)
        click.echo("Pipeline Summary")
        click.echo("=" * 60)

        if trim_stats:
            click.echo(f"\nAdapter Trimming: {trim_stats.get('percent_trimmed', 0):.2f}% trimmed")

        click.echo(f"Quality Trimming: {qtrim_stats.get('percent_trimmed', 0):.2f}% trimmed")
        click.echo(f"Filtering: {filter_stats.get('percent_passed', 0):.2f}% passed")

        # Quality scores
        metrics = QualityMetrics()
        initial_score = metrics.calculate_composite_score(initial_qc)
        final_score = metrics.calculate_composite_score(final_qc)

        click.echo(f"\nInitial Quality Score: {initial_score}/100")
        click.echo(f"Final Quality Score: {final_score}/100")
        click.echo(f"Improvement: {final_score - initial_score:+.2f}")

        click.echo(f"\nFinal output: {filtered_file}")
        click.echo(f"QC reports: {output_dir}")

    # Cleanup intermediate files if requested
    if not keep_intermediate and adapter:
        Path(trimmed_file).unlink(missing_ok=True)
        Path(qtrimmed_file).unlink(missing_ok=True)


@cli.command()
@click.option(
    "--sample-sheet",
    "-s",
    required=True,
    type=click.Path(exists=True),
    help="Sample sheet CSV file",
)
@click.option(
    "--output",
    "-o",
    "output_dir",
    required=True,
    type=click.Path(),
    help="Output directory",
)
@click.option("--threads", "-t", default=4, help="Number of parallel workers")
@click.option("--quiet", "-q", is_flag=True, help="Suppress output")
def batch(sample_sheet: str, output_dir: str, threads: int, quiet: bool) -> None:
    """Process multiple samples from a sample sheet."""
    if not quiet:
        click.echo(f"Processing samples from {sample_sheet}...")

    # Parse sample sheet
    samples = parse_sample_sheet(sample_sheet)

    if not quiet:
        click.echo(f"Found {len(samples)} samples to process")

    def process_sample(sample: dict[str, str]) -> dict[str, str]:
        """Process a single sample."""
        sample_id = sample["sample_id"]
        output_path = Path(output_dir) / sample_id

        # Build config
        config = Config(
            min_quality=int(sample.get("min_quality", 20)),
            min_length=int(sample.get("min_length", 50)),
            adapter_sequence=sample.get("adapter"),
            output_dir=output_path,
        )

        # Run QC
        runner = FastQCRunner(config)
        qc_dir = output_path / "qc"
        runner.run(sample["fastq_r1"], str(qc_dir), quiet=True)

        return {"sample_id": sample_id, "status": "completed"}

    # Process in parallel
    results = process_samples_parallel(
        samples,
        process_sample,
        n_workers=threads,
        show_progress=not quiet,
        description="Processing samples",
    )

    if not quiet:
        completed = sum(1 for r in results if r and r.get("status") == "completed")
        click.echo(f"\nProcessed {completed}/{len(samples)} samples successfully")
        click.echo(f"Results saved to: {output_dir}")


@cli.command(name="validate-pairs")
@click.option(
    "--r1",
    required=True,
    type=click.Path(exists=True),
    help="R1 FASTQ file",
)
@click.option(
    "--r2",
    required=True,
    type=click.Path(exists=True),
    help="R2 FASTQ file",
)
@click.option(
    "--max-pairs",
    type=int,
    help="Maximum number of pairs to check (default: all)",
)
@click.option("--no-order-check", is_flag=True, help="Skip read order validation")
@click.option("--quiet", "-q", is_flag=True, help="Suppress output")
def validate_pairs(
    r1: str,
    r2: str,
    max_pairs: int | None,
    no_order_check: bool,
    quiet: bool,
) -> None:
    """
    Validate paired-end FASTQ files (R1/R2 concordance).

    Checks for:
      - Matching read counts
      - Correct read ID pairing
      - File integrity

    Examples:
      bioseqflow validate-pairs --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz
    """
    if not quiet:
        echo_header("Paired-End Validation")
        echo_info(f"R1 file: {r1}")
        echo_info(f"R2 file: {r2}")

    result = validate_paired_files(
        r1,
        r2,
        check_order=not no_order_check,
        max_pairs_to_check=max_pairs,
    )

    if not quiet:
        click.echo()
        if result.is_valid:
            echo_success(f"Validation passed! {result.total_pairs:,} pairs validated")
        else:
            echo_error("Validation failed!")
            for error in result.errors:
                click.echo(click.style(f"  - {error}", fg="red"))
            click.get_current_context().exit(1)


@cli.command(name="insert-size")
@click.option(
    "--r1",
    required=True,
    type=click.Path(exists=True),
    help="R1 FASTQ file",
)
@click.option(
    "--r2",
    required=True,
    type=click.Path(exists=True),
    help="R2 FASTQ file",
)
@click.option(
    "--max-reads",
    default=10000,
    type=int,
    help="Maximum reads to analyze (default: 10000)",
)
def insert_size(r1: str, r2: str, max_reads: int) -> None:
    """
    Calculate insert size distribution for paired-end data.

    Examples:
      bioseqflow insert-size --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz
    """
    echo_header("Insert Size Analysis")
    echo_info(f"Analyzing {max_reads:,} read pairs...")

    stats = calculate_insert_size_distribution(r1, r2, max_reads=max_reads)

    click.echo()
    click.echo(click.style("Insert Size Statistics:", fg="cyan", bold=True))
    click.echo(f"  Pairs analyzed: {stats['pairs_analyzed']:,}")
    click.echo(
        "  Mean: " + click.style(f"{stats['mean_insert_size']:.1f} bp", fg="yellow", bold=True)
    )
    click.echo(f"  Median: {stats['median_insert_size']:.1f} bp")
    click.echo(f"  Std dev: {stats['std_insert_size']:.1f} bp")
    click.echo(f"  Range: {stats['min_insert_size']}-{stats['max_insert_size']} bp")


@cli.command(name="check-orientation")
@click.option(
    "--r1",
    required=True,
    type=click.Path(exists=True),
    help="R1 FASTQ file",
)
@click.option(
    "--r2",
    required=True,
    type=click.Path(exists=True),
    help="R2 FASTQ file",
)
@click.option(
    "--max-reads",
    default=1000,
    type=int,
    help="Maximum reads to check (default: 1000)",
)
def check_orientation_cmd(r1: str, r2: str, max_reads: int) -> None:
    """
    Check read orientation and length consistency.

    Examples:
      bioseqflow check-orientation --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz
    """
    echo_header("Read Orientation Check")

    stats = check_read_orientation(r1, r2, max_reads=max_reads)

    click.echo()
    click.echo(click.style("Orientation Statistics:", fg="cyan", bold=True))
    click.echo(f"  Pairs checked: {stats['pairs_checked']:,}")
    click.echo(f"  Mean R1 length: {stats['mean_r1_length']:.1f} bp")
    click.echo(f"  Mean R2 length: {stats['mean_r2_length']:.1f} bp")
    click.echo(f"  Length difference: {stats['mean_length_diff']:.1f} bp")

    click.echo()
    if stats["length_consistent"]:
        echo_success("Read lengths are consistent (diff < 5 bp)")
    else:
        echo_warning(
            f"Read lengths differ significantly (mean diff: {stats['mean_length_diff']:.1f} bp)"
        )


@cli.command(name="parse-qc")
@click.option(
    "--fastqc",
    type=click.Path(exists=True),
    help="FastQC data file (fastqc_data.txt)",
)
@click.option(
    "--multiqc",
    type=click.Path(exists=True),
    help="MultiQC data directory",
)
def parse_qc(fastqc: str | None, multiqc: str | None) -> None:
    """
    Parse and display QC results from FastQC or MultiQC.

    Examples:
      # Parse FastQC output
      bioseqflow parse-qc --fastqc sample_fastqc/fastqc_data.txt

      # Parse MultiQC output
      bioseqflow parse-qc --multiqc multiqc_data/
    """
    if not fastqc and not multiqc:
        echo_error("Please provide either --fastqc or --multiqc")
        click.get_current_context().exit(1)

    if fastqc:
        echo_header("FastQC Results")
        parser = FastQCParser()
        results = parser.parse(fastqc)

        # Basic statistics
        if "Basic Statistics" in results["modules"]:
            stats = results["modules"]["Basic Statistics"]
            click.echo(click.style("\nBasic Statistics:", fg="cyan", bold=True))
            for key, value in stats.items():
                click.echo(f"  {key}: {value}")

        # Module status
        click.echo(click.style("\nModule Status:", fg="cyan", bold=True))
        for module, status in results["summary"].items():
            if status == "pass":
                click.echo(click.style(f"  [PASS] {module}", fg="green"))
            elif status == "warn":
                click.echo(click.style(f"  [WARN] {module}", fg="yellow"))
            else:
                click.echo(click.style(f"  [FAIL] {module}", fg="red"))

    if multiqc:
        echo_header("MultiQC Results")
        multiqc_parser = MultiQCParser()
        results = multiqc_parser.parse(multiqc)

        if results["general_stats"]:
            click.echo(click.style("\nSample Statistics:", fg="cyan", bold=True))
            for sample, stats in list(results["general_stats"].items())[:10]:
                click.echo(click.style(f"\n  {sample}", fg="yellow", bold=True))
                for key, value in list(stats.items())[:5]:
                    click.echo(f"    {key}: {value}")

            if len(results["general_stats"]) > 10:
                echo_info(f"... and {len(results['general_stats']) - 10} more samples")


if __name__ == "__main__":
    cli()
