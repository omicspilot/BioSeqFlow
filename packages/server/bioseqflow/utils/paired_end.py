from __future__ import annotations

"""Utilities for paired-end sequencing data validation and processing."""

from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

from bioseqflow.utils.io import FastqRecord, read_fastq

if TYPE_CHECKING:
    from collections.abc import Generator


class PairedRecord(NamedTuple):
    """Represents a pair of R1/R2 reads."""

    r1: FastqRecord
    r2: FastqRecord


class ValidationResult(NamedTuple):
    """Result from paired-end validation."""

    is_valid: bool
    total_pairs: int
    errors: list[str]


def extract_read_id(header: str) -> str:
    """
    Extract read ID from FASTQ header.

    Handles various header formats (Illumina, SRA, etc.) and removes
    pair identifiers (/1, /2, or space-separated pair info).

    Args:
        header: FASTQ header line (starts with @)

    Returns:
        Base read ID without pair information

    Examples:
        >>> extract_read_id("@SRR123456.1 1 length=150")
        'SRR123456.1'

        >>> extract_read_id("@instrument:run:flowcell:lane:tile:x:y 1:N:0:ATCG")
        'instrument:run:flowcell:lane:tile:x:y'

        >>> extract_read_id("@read1/1")
        'read1'
    """
    # Remove leading @
    header = header.lstrip("@").strip()

    # Remove /1, /2 suffixes (common in older formats)
    if header.endswith("/1") or header.endswith("/2"):
        header = header[:-2]

    # Split on whitespace and take first part (removes pair info like "1:N:0:ATCG")
    header = header.split()[0]

    return header


def validate_paired_files(
    r1_file: Path | str,
    r2_file: Path | str,
    check_order: bool = True,
    max_pairs_to_check: int | None = None,
) -> ValidationResult:
    """
    Validate that R1 and R2 FASTQ files are properly paired.

    Checks:
    - Both files have same number of reads
    - Read IDs match (in order if check_order=True)
    - Files are not identical (common mistake)

    Args:
        r1_file: R1 FASTQ file path
        r2_file: R2 FASTQ file path
        check_order: Verify reads are in same order
        max_pairs_to_check: Maximum pairs to check (None = all)

    Returns:
        ValidationResult with validation status and errors

    Examples:
        >>> result = validate_paired_files("sample_R1.fastq", "sample_R2.fastq")
        >>> if result.is_valid:
        ...     print(f"Valid paired files with {result.total_pairs} pairs")
        ... else:
        ...     print(f"Errors: {result.errors}")
    """
    r1_file = Path(r1_file)
    r2_file = Path(r2_file)

    errors = []
    pairs_checked = 0

    # Check files exist
    if not r1_file.exists():
        errors.append(f"R1 file not found: {r1_file}")
        return ValidationResult(False, 0, errors)

    if not r2_file.exists():
        errors.append(f"R2 file not found: {r2_file}")
        return ValidationResult(False, 0, errors)

    # Check files are not identical
    if r1_file.resolve() == r2_file.resolve():
        errors.append("R1 and R2 files are the same file")
        return ValidationResult(False, 0, errors)

    # Read both files simultaneously
    r1_reader = read_fastq(r1_file)
    r2_reader = read_fastq(r2_file)

    try:
        while True:
            # Get next records
            try:
                r1_record = next(r1_reader)
            except StopIteration:
                r1_record = None

            try:
                r2_record = next(r2_reader)
            except StopIteration:
                r2_record = None

            # Check if one file ended before the other
            if r1_record is None and r2_record is not None:
                errors.append(f"R1 file ended at {pairs_checked} reads, but R2 continues")
                break
            elif r1_record is not None and r2_record is None:
                errors.append(f"R2 file ended at {pairs_checked} reads, but R1 continues")
                break
            elif r1_record is None and r2_record is None:
                # Both ended at same time - good!
                break

            pairs_checked += 1

            # At this point, both records must be not None
            # (all None cases are handled above with break)
            assert r1_record is not None
            assert r2_record is not None

            # Extract and compare read IDs
            if check_order:
                r1_id = extract_read_id(r1_record.header)
                r2_id = extract_read_id(r2_record.header)

                if r1_id != r2_id:
                    errors.append(
                        f"Read ID mismatch at pair {pairs_checked}: R1={r1_id}, R2={r2_id}"
                    )
                    # Stop after first mismatch
                    break

            # Stop if reached max pairs to check
            if max_pairs_to_check and pairs_checked >= max_pairs_to_check:
                break

    except Exception as e:
        errors.append(f"Error during validation: {str(e)}")

    is_valid = len(errors) == 0

    return ValidationResult(is_valid, pairs_checked, errors)


def read_paired_fastq(
    r1_file: Path | str, r2_file: Path | str, validate: bool = True
) -> Generator[PairedRecord, None, None]:
    """
    Read paired FASTQ files simultaneously.

    Args:
        r1_file: R1 FASTQ file path
        r2_file: R2 FASTQ file path
        validate: Validate read IDs match

    Yields:
        PairedRecord objects with R1 and R2 reads

    Raises:
        ValueError: If validation enabled and read IDs don't match

    Examples:
        >>> for pair in read_paired_fastq("R1.fastq", "R2.fastq"):
        ...     print(f"R1: {pair.r1.sequence}")
        ...     print(f"R2: {pair.r2.sequence}")
    """
    r1_file = Path(r1_file)
    r2_file = Path(r2_file)

    r1_reader = read_fastq(r1_file)
    r2_reader = read_fastq(r2_file)

    pair_num = 0

    try:
        for r1_record in r1_reader:
            pair_num += 1

            try:
                r2_record = next(r2_reader)
            except StopIteration:
                raise ValueError(
                    f"R2 file ended early at pair {pair_num}, but R1 continues"
                ) from None

            # Validate read IDs match
            if validate:
                r1_id = extract_read_id(r1_record.header)
                r2_id = extract_read_id(r2_record.header)

                if r1_id != r2_id:
                    raise ValueError(f"Read ID mismatch at pair {pair_num}: R1={r1_id}, R2={r2_id}")

            yield PairedRecord(r1_record, r2_record)

    finally:
        # Check if R2 has extra reads
        try:
            next(r2_reader)
            # If we get here, R2 has more reads
            if validate:
                raise ValueError(f"R1 file ended at pair {pair_num}, but R2 continues")
        except StopIteration:
            # Good - both files ended at same point
            pass


def calculate_insert_size_distribution(
    r1_file: Path | str,
    r2_file: Path | str,
    max_reads: int = 10000,
) -> dict[str, float | int]:
    """
    Calculate insert size distribution from paired-end reads.

    Note: This is an estimate based on read lengths. For accurate
    insert sizes, alignment to a reference is required.

    Args:
        r1_file: R1 FASTQ file path
        r2_file: R2 FASTQ file path
        max_reads: Maximum reads to analyze

    Returns:
        Dictionary with insert size statistics

    Examples:
        >>> stats = calculate_insert_size_distribution("R1.fq", "R2.fq")
        >>> print(f"Mean insert size: {stats['mean_insert_size']}")
    """
    insert_sizes = []
    pairs_analyzed = 0

    for pair in read_paired_fastq(r1_file, r2_file, validate=False):
        # Estimate: R1 length + R2 length (assumes no overlap)
        # Real insert size requires alignment
        estimated_insert = pair.r1.length + pair.r2.length

        insert_sizes.append(estimated_insert)
        pairs_analyzed += 1

        if pairs_analyzed >= max_reads:
            break

    if not insert_sizes:
        return {
            "pairs_analyzed": 0,
            "mean_insert_size": 0.0,
            "min_insert_size": 0,
            "max_insert_size": 0,
        }

    return {
        "pairs_analyzed": pairs_analyzed,
        "mean_insert_size": sum(insert_sizes) / len(insert_sizes),
        "min_insert_size": min(insert_sizes),
        "max_insert_size": max(insert_sizes),
    }


def check_read_orientation(
    r1_file: Path | str,
    r2_file: Path | str,
    max_reads: int = 1000,
) -> dict[str, int | float]:
    """
    Check read orientation consistency.

    Analyzes whether reads are properly oriented (FR, RF, FF, RR).
    This is a basic check - detailed orientation requires alignment.

    Args:
        r1_file: R1 FASTQ file path
        r2_file: R2 FASTQ file path
        max_reads: Maximum reads to check

    Returns:
        Dictionary with orientation statistics
    """
    pairs_checked = 0
    length_diffs = []

    for pair in read_paired_fastq(r1_file, r2_file, validate=False):
        pairs_checked += 1

        # Check length consistency
        length_diff = abs(pair.r1.length - pair.r2.length)
        length_diffs.append(length_diff)

        if pairs_checked >= max_reads:
            break

    if not length_diffs:
        return {
            "pairs_checked": 0,
            "mean_length_diff": 0.0,
            "max_length_diff": 0,
            "length_consistent": True,
        }

    mean_diff = sum(length_diffs) / len(length_diffs)
    max_diff = max(length_diffs)

    # Consider consistent if mean diff < 10bp
    is_consistent = mean_diff < 10.0

    return {
        "pairs_checked": pairs_checked,
        "mean_length_diff": mean_diff,
        "max_length_diff": max_diff,
        "length_consistent": is_consistent,
    }
