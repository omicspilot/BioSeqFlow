#!/usr/bin/env python3
"""Example: Parsing QC tool outputs (FastQC and MultiQC).

This example demonstrates:
1. Parsing FastQC data files
2. Extracting quality metrics
3. Parsing MultiQC aggregated results
4. Generating QC summaries
"""

from pathlib import Path

from bioseqflow.qc.parsers import FastQCParser, MultiQCParser


def create_example_fastqc_data(output_dir: Path) -> Path:
    """Create an example FastQC data file."""
    output_dir.mkdir(exist_ok=True, parents=True)
    fastqc_file = output_dir / "fastqc_data.txt"

    content = """##FastQC	0.11.9
>>Basic Statistics	pass
#Measure	Value
Filename	sample_R1.fastq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	5000000
Sequences flagged as poor quality	0
Sequence length	150
%GC	45
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	37.5	38	37	38	36	39
2	37.3	38	37	38	35	39
3	37.1	37	37	38	35	38
>>END_MODULE
>>Per sequence quality scores	warn
#Quality	Count
30	150000
35	1850000
40	3000000
>>END_MODULE
>>Per base sequence content	fail
#Base	G	A	T	C
1	25.2	24.8	24.5	25.5
2	25.1	24.9	24.6	25.4
>>END_MODULE
>>Sequence Length Distribution	pass
#Length	Count
150	5000000
>>END_MODULE
>>Adapter Content	pass
#Position	Illumina Universal Adapter	Illumina Small RNA 3' Adapter	Nextera Transposase Sequence
1	0.0	0.0	0.0
2	0.0	0.0	0.0
>>END_MODULE
"""
    fastqc_file.write_text(content)
    return fastqc_file


def create_example_multiqc_data(output_dir: Path) -> Path:
    """Create example MultiQC data directory."""
    multiqc_dir = output_dir / "multiqc_data"
    multiqc_dir.mkdir(exist_ok=True, parents=True)

    # General stats file
    general_stats = multiqc_dir / "multiqc_general_stats.txt"
    general_stats.write_text("""Sample	Total Sequences	%GC	avg_quality	percent_duplicates
sample1_R1	5000000	45.2	37.5	12.3
sample1_R2	5000000	44.8	37.2	12.1
sample2_R1	7500000	42.1	36.8	15.2
sample2_R2	7500000	41.9	36.5	15.4
""")

    # FastQC data file
    fastqc_data = multiqc_dir / "multiqc_fastqc.txt"
    fastqc_data.write_text("""Sample	total_sequences	percent_gc	avg_sequence_length	percent_duplicates	percent_fails
sample1_R1	5000000	45.2	150	12.3	0.5
sample1_R2	5000000	44.8	150	12.1	0.4
sample2_R1	7500000	42.1	150	15.2	0.8
sample2_R2	7500000	41.9	150	15.4	0.7
""")

    return multiqc_dir


def example_parse_fastqc():
    """Demonstrate parsing FastQC output."""
    print("=" * 60)
    print("Parsing FastQC Output")
    print("=" * 60)

    # Create example data
    example_dir = Path("/tmp/qc_parsing_example")
    fastqc_file = create_example_fastqc_data(example_dir)

    # Parse FastQC data
    parser = FastQCParser()
    results = parser.parse(fastqc_file)

    # Display basic statistics
    print("\n--- Basic Statistics ---")
    basic_stats = results["modules"]["Basic Statistics"]
    for key, value in basic_stats.items():
        print(f"  {key}: {value}")

    # Display module status summary
    print("\n--- Module Status Summary ---")
    for module, status in results["summary"].items():
        status_icon = {
            "pass": "✓",
            "warn": "⚠",
            "fail": "✗"
        }.get(status, "?")
        print(f"  {status_icon} {module}: {status.upper()}")

    # Analyze quality
    print("\n--- Quality Analysis ---")
    fail_count = sum(1 for status in results["summary"].values() if status == "fail")
    warn_count = sum(1 for status in results["summary"].values() if status == "warn")
    pass_count = sum(1 for status in results["summary"].values() if status == "pass")

    print(f"  Passed: {pass_count} modules")
    print(f"  Warnings: {warn_count} modules")
    print(f"  Failed: {fail_count} modules")

    # Overall assessment
    if fail_count == 0 and warn_count == 0:
        print("\n  Overall: ✓ EXCELLENT quality")
    elif fail_count == 0:
        print("\n  Overall: ⚠ GOOD quality (some warnings)")
    else:
        print("\n  Overall: ✗ NEEDS ATTENTION (failures detected)")


def example_parse_multiqc():
    """Demonstrate parsing MultiQC output."""
    print("\n" + "=" * 60)
    print("Parsing MultiQC Output")
    print("=" * 60)

    # Create example data
    example_dir = Path("/tmp/qc_parsing_example")
    multiqc_dir = create_example_multiqc_data(example_dir)

    # Parse MultiQC data
    parser = MultiQCParser()
    results = parser.parse(multiqc_dir)

    # Display general statistics
    print("\n--- General Statistics (All Samples) ---")
    print(f"{'Sample':<15} {'Total Seqs':>12} {'%GC':>6} {'Avg Qual':>9} {'% Dup':>7}")
    print("-" * 60)

    for sample, stats in results["general_stats"].items():
        print(
            f"{sample:<15} "
            f"{stats['Total Sequences']:>12} "
            f"{stats['%GC']:>6} "
            f"{stats['avg_quality']:>9} "
            f"{stats['percent_duplicates']:>7}"
        )

    # Display FastQC-specific data
    print("\n--- FastQC Data (All Samples) ---")
    if results["fastqc_data"]:
        print(f"{'Sample':<15} {'Sequences':>12} {'GC%':>6} {'Length':>7} {'Dup%':>6} {'Fail%':>6}")
        print("-" * 60)

        for sample, data in results["fastqc_data"].items():
            print(
                f"{sample:<15} "
                f"{data['total_sequences']:>12} "
                f"{data['percent_gc']:>6} "
                f"{data['avg_sequence_length']:>7} "
                f"{data['percent_duplicates']:>6} "
                f"{data['percent_fails']:>6}"
            )

    # Summary statistics
    print("\n--- Summary Across All Samples ---")
    if results["general_stats"]:
        samples = list(results["general_stats"].keys())
        total_seqs = sum(
            int(stats["Total Sequences"])
            for stats in results["general_stats"].values()
        )
        avg_gc = sum(
            float(stats["%GC"])
            for stats in results["general_stats"].values()
        ) / len(samples)

        print(f"  Total samples: {len(samples)}")
        print(f"  Total sequences: {total_seqs:,}")
        print(f"  Average GC content: {avg_gc:.1f}%")


def example_compare_samples():
    """Demonstrate comparing samples from MultiQC output."""
    print("\n" + "=" * 60)
    print("Sample Comparison")
    print("=" * 60)

    # Create example data
    example_dir = Path("/tmp/qc_parsing_example")
    multiqc_dir = create_example_multiqc_data(example_dir)

    # Parse data
    parser = MultiQCParser()
    results = parser.parse(multiqc_dir)

    # Compare R1 vs R2 for each sample
    print("\n--- R1 vs R2 Comparison ---")

    # Group by sample
    samples = {}
    for sample_name in results["fastqc_data"].keys():
        base_name = sample_name.rsplit("_", 1)[0]  # Remove _R1/_R2
        if base_name not in samples:
            samples[base_name] = {}
        if "_R1" in sample_name:
            samples[base_name]["R1"] = results["fastqc_data"][sample_name]
        elif "_R2" in sample_name:
            samples[base_name]["R2"] = results["fastqc_data"][sample_name]

    # Compare
    for sample_name, reads in samples.items():
        if "R1" in reads and "R2" in reads:
            print(f"\n{sample_name}:")
            r1_gc = float(reads["R1"]["percent_gc"])
            r2_gc = float(reads["R2"]["percent_gc"])
            gc_diff = abs(r1_gc - r2_gc)

            print(f"  GC content: R1={r1_gc:.1f}%, R2={r2_gc:.1f}% (diff={gc_diff:.1f}%)")

            if gc_diff < 2.0:
                print("  ✓ GC content is consistent between R1 and R2")
            else:
                print("  ⚠ GC content differs significantly between R1 and R2")


def example_filter_low_quality_samples():
    """Demonstrate filtering samples based on QC metrics."""
    print("\n" + "=" * 60)
    print("Quality Filtering")
    print("=" * 60)

    # Create example data
    example_dir = Path("/tmp/qc_parsing_example")
    multiqc_dir = create_example_multiqc_data(example_dir)

    # Parse data
    parser = MultiQCParser()
    results = parser.parse(multiqc_dir)

    # Define quality thresholds
    MIN_AVG_QUALITY = 35.0
    MAX_DUPLICATION = 20.0
    MAX_FAILS = 1.0

    print("\nQuality Thresholds:")
    print(f"  Minimum average quality: {MIN_AVG_QUALITY}")
    print(f"  Maximum duplication rate: {MAX_DUPLICATION}%")
    print(f"  Maximum failure rate: {MAX_FAILS}%")

    print("\n--- Sample Quality Assessment ---")

    passed_samples = []
    failed_samples = []

    for sample, stats in results["general_stats"].items():
        avg_qual = float(stats["avg_quality"])
        dup_rate = float(stats["percent_duplicates"])

        # Get fail rate from FastQC data if available
        fail_rate = 0.0
        if sample in results["fastqc_data"]:
            fail_rate = float(results["fastqc_data"][sample]["percent_fails"])

        # Check quality
        passes = True
        issues = []

        if avg_qual < MIN_AVG_QUALITY:
            passes = False
            issues.append(f"low quality ({avg_qual:.1f})")

        if dup_rate > MAX_DUPLICATION:
            passes = False
            issues.append(f"high duplication ({dup_rate:.1f}%)")

        if fail_rate > MAX_FAILS:
            passes = False
            issues.append(f"high failure rate ({fail_rate:.1f}%)")

        if passes:
            passed_samples.append(sample)
            print(f"  ✓ {sample}: PASS")
        else:
            failed_samples.append(sample)
            print(f"  ✗ {sample}: FAIL - {', '.join(issues)}")

    print(f"\nSummary:")
    print(f"  Passed: {len(passed_samples)}/{len(results['general_stats'])} samples")
    print(f"  Failed: {len(failed_samples)}/{len(results['general_stats'])} samples")

    if failed_samples:
        print(f"\n  Failed samples: {', '.join(failed_samples)}")
        print("  → These samples may need resequencing or should be excluded from analysis")


def main():
    """Run all QC parsing examples."""
    example_parse_fastqc()
    example_parse_multiqc()
    example_compare_samples()
    example_filter_low_quality_samples()

    print("\n" + "=" * 60)
    print("QC Parsing Best Practices:")
    print("=" * 60)
    print("1. Always check module status (pass/warn/fail) from FastQC")
    print("2. Compare R1 vs R2 for paired-end data consistency")
    print("3. Set quality thresholds appropriate for your experiment")
    print("4. Use MultiQC to aggregate results across many samples")
    print("5. Flag samples with high duplication or failure rates")
    print("6. Document QC thresholds in your analysis pipeline")
    print("=" * 60)


if __name__ == "__main__":
    main()
