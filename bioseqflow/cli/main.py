from __future__ import annotations

"""Command-line interface for BioSeqFlow."""

from pathlib import Path
from typing import Optional

import click

from bioseqflow import __version__
from bioseqflow.core.config import Config
from bioseqflow.preprocessing.deduplication import DuplicateRemover
from bioseqflow.preprocessing.filtering import QualityFilter
from bioseqflow.preprocessing.trimming import AdapterTrimmer, QualityTrimmer
from bioseqflow.qc.fastqc import FastQCRunner
from bioseqflow.qc.metrics import QualityMetrics
from bioseqflow.utils.io import parse_sample_sheet
from bioseqflow.utils.parallel import process_samples_parallel


@click.group()
@click.version_option(version=__version__)
def cli() -> None:
    """BioSeqFlow: Sequencing Quality Control and Preprocessing Platform."""
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
    """Run quality control analysis."""
    click.echo(f"Running quality control on {input_file}...")

    runner = FastQCRunner()
    results = runner.run(input_file, output_dir, threads=threads, quiet=quiet)

    if not quiet:
        click.echo("\n=== QC Results ===")
        if "basic_statistics" in results.get("modules", {}):
            stats = results["modules"]["basic_statistics"]
            click.echo(f"Total Sequences: {stats.get('Total Sequences', 'N/A')}")
            click.echo(f"Sequence Length: {stats.get('Sequence length', 'N/A')}")
            click.echo(f"%GC: {stats.get('%GC', 'N/A')}")

        # Calculate composite score
        metrics = QualityMetrics()
        score = metrics.calculate_composite_score(results)
        click.echo(f"\nComposite Quality Score: {score}/100")

        if score >= 90:
            click.echo(click.style("Quality: Excellent", fg="green"))
        elif score >= 70:
            click.echo(click.style("Quality: Good", fg="green"))
        elif score >= 50:
            click.echo(click.style("Quality: Needs Review", fg="yellow"))
        else:
            click.echo(click.style("Quality: Requires Attention", fg="red"))

    click.echo(f"\nResults saved to: {output_dir}")


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
@click.option("--adapter", "-a", required=True, help="Adapter sequence")
@click.option("--quiet", "-q", is_flag=True, help="Suppress output")
def trim(input_file: str, output_file: str, adapter: str, quiet: bool) -> None:
    """Trim adapter sequences."""
    if not quiet:
        click.echo(f"Trimming adapters from {input_file}...")

    trimmer = AdapterTrimmer()
    stats = trimmer.trim(input_file, output_file, adapter)

    if not quiet:
        click.echo("\n=== Trimming Results ===")
        click.echo(f"Total Reads: {stats.get('total_reads', 'N/A')}")
        click.echo(f"Trimmed Reads: {stats.get('trimmed_reads', 'N/A')}")
        click.echo(f"Percent Trimmed: {stats.get('percent_trimmed', 0):.2f}%")
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
    adapter: Optional[str],
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


if __name__ == "__main__":
    cli()
