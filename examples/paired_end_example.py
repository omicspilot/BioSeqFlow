#!/usr/bin/env python3
"""Example: Working with paired-end sequencing data.

This example demonstrates:
1. Validating paired-end FASTQ files
2. Reading paired reads simultaneously
3. Calculating insert size distribution
4. Checking read orientation consistency
"""

from pathlib import Path

from bioseqflow.utils.io import FastqRecord, write_fastq
from bioseqflow.utils.paired_end import (
    calculate_insert_size_distribution,
    check_read_orientation,
    extract_read_id,
    read_paired_fastq,
    validate_paired_files,
)


def create_example_paired_files(output_dir: Path) -> tuple[Path, Path]:
    """Create example paired-end FASTQ files for demonstration."""
    output_dir.mkdir(exist_ok=True)

    r1_file = output_dir / "example_R1.fastq"
    r2_file = output_dir / "example_R2.fastq"

    # Create R1 reads (forward)
    r1_records = [
        FastqRecord(
            "@SEQ_1 1:N:0:ATCG\n",
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
        FastqRecord(
            "@SEQ_2 1:N:0:ATCG\n",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
        FastqRecord(
            "@SEQ_3 1:N:0:ATCG\n",
            "TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
    ]

    # Create R2 reads (reverse) - same read IDs
    r2_records = [
        FastqRecord(
            "@SEQ_1 2:N:0:ATCG\n",
            "CGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACG\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
        FastqRecord(
            "@SEQ_2 2:N:0:ATCG\n",
            "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
        FastqRecord(
            "@SEQ_3 2:N:0:ATCG\n",
            "AATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATT\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
    ]

    write_fastq(r1_records, r1_file)
    write_fastq(r2_records, r2_file)

    return r1_file, r2_file


def example_read_id_extraction():
    """Demonstrate read ID extraction from various FASTQ header formats."""
    print("=" * 60)
    print("Read ID Extraction Examples")
    print("=" * 60)

    headers = [
        "@SEQ_1",  # Simple format
        "@instrument:run:flowcell:lane:tile:x:y/1",  # Old Illumina format
        "@instrument:run:flowcell:lane:tile:x:y 1:N:0:ATCG",  # New Illumina
        "@SRR123456.1 1 length=150",  # SRA format
    ]

    for header in headers:
        read_id = extract_read_id(header)
        print(f"Header:  {header}")
        print(f"Read ID: {read_id}\n")


def example_validate_paired_files(r1_file: Path, r2_file: Path):
    """Demonstrate paired-end file validation."""
    print("\n" + "=" * 60)
    print("Paired-End File Validation")
    print("=" * 60)

    # Validate files
    result = validate_paired_files(r1_file, r2_file, check_order=True)

    print(f"\nR1 file: {r1_file.name}")
    print(f"R2 file: {r2_file.name}")
    print(f"Valid: {result.is_valid}")
    print(f"Total pairs checked: {result.total_pairs}")

    if not result.is_valid:
        print("\nValidation errors:")
        for error in result.errors:
            print(f"  - {error}")
    else:
        print("✓ All validation checks passed!")

    # Example: Validate with limited pairs
    print("\n" + "-" * 60)
    print("Validating first 1000 pairs only...")
    quick_result = validate_paired_files(
        r1_file, r2_file, max_pairs_to_check=1000
    )
    print(f"Validated {quick_result.total_pairs} pairs")


def example_read_paired_fastq(r1_file: Path, r2_file: Path):
    """Demonstrate reading paired FASTQ files."""
    print("\n" + "=" * 60)
    print("Reading Paired FASTQ Files")
    print("=" * 60)

    pair_count = 0
    for pair in read_paired_fastq(r1_file, r2_file, validate=True):
        pair_count += 1
        if pair_count <= 2:  # Show first 2 pairs
            print(f"\nPair {pair_count}:")
            print(f"  R1 header: {pair.r1.header.strip()}")
            print(f"  R1 length: {pair.r1.length} bp")
            print(f"  R2 header: {pair.r2.header.strip()}")
            print(f"  R2 length: {pair.r2.length} bp")

    print(f"\nTotal pairs processed: {pair_count}")


def example_insert_size_distribution(r1_file: Path, r2_file: Path):
    """Demonstrate insert size calculation."""
    print("\n" + "=" * 60)
    print("Insert Size Distribution")
    print("=" * 60)

    stats = calculate_insert_size_distribution(r1_file, r2_file, max_reads=1000)

    print(f"\nPairs analyzed: {stats['pairs_analyzed']}")
    print(f"Mean insert size: {stats['mean_insert_size']:.1f} bp")
    print(f"Median insert size: {stats['median_insert_size']:.1f} bp")
    print(f"Std dev: {stats['std_insert_size']:.1f} bp")
    print(f"Min insert size: {stats['min_insert_size']} bp")
    print(f"Max insert size: {stats['max_insert_size']} bp")


def example_read_orientation(r1_file: Path, r2_file: Path):
    """Demonstrate read orientation checking."""
    print("\n" + "=" * 60)
    print("Read Orientation Check")
    print("=" * 60)

    stats = check_read_orientation(r1_file, r2_file, max_reads=1000)

    print(f"\nPairs checked: {stats['pairs_checked']}")
    print(f"Mean R1 length: {stats['mean_r1_length']:.1f} bp")
    print(f"Mean R2 length: {stats['mean_r2_length']:.1f} bp")
    print(f"Mean length difference: {stats['mean_length_diff']:.1f} bp")
    print(f"Max length difference: {stats['max_length_diff']} bp")

    if stats['length_consistent']:
        print("\n✓ Read lengths are consistent (mean diff < 5 bp)")
    else:
        print("\n⚠ Read lengths show inconsistency (mean diff >= 5 bp)")


def main():
    """Run all paired-end examples."""
    # Setup
    example_dir = Path("/tmp/paired_end_example")
    print("Creating example paired-end FASTQ files...")
    r1_file, r2_file = create_example_paired_files(example_dir)

    # Run examples
    example_read_id_extraction()
    example_validate_paired_files(r1_file, r2_file)
    example_read_paired_fastq(r1_file, r2_file)
    example_insert_size_distribution(r1_file, r2_file)
    example_read_orientation(r1_file, r2_file)

    print("\n" + "=" * 60)
    print(f"Example files saved to: {example_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
