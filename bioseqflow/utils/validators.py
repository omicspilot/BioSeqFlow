from __future__ import annotations

"""Validation utilities for input checking."""

from pathlib import Path


def validate_file_path(
    file_path: Path | str, must_exist: bool = True, extensions: list[str] | None = None
) -> Path:
    """
    Validate file path.

    Args:
        file_path: Path to validate
        must_exist: Whether file must exist
        extensions: Allowed file extensions (e.g., ['.fastq', '.fq', '.fastq.gz'])

    Returns:
        Validated Path object

    Raises:
        FileNotFoundError: If file doesn't exist and must_exist=True
        ValueError: If file extension is not allowed
    """
    file_path = Path(file_path)

    if must_exist and not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    if extensions:
        # Handle double extensions like .fastq.gz
        full_suffix = "".join(file_path.suffixes)
        single_suffix = file_path.suffix

        if full_suffix not in extensions and single_suffix not in extensions:
            raise ValueError(
                f"Invalid file extension. Expected one of {extensions}, got {full_suffix}"
            )

    return file_path


def validate_quality_score(quality: int, encoding: str = "phred33") -> None:
    """
    Validate quality score.

    Args:
        quality: Quality score to validate
        encoding: Quality encoding ('phred33' or 'phred64')

    Raises:
        ValueError: If quality score is out of range
    """
    if encoding == "phred33":
        if not 0 <= quality <= 93:
            raise ValueError(f"Phred+33 quality must be 0-93, got {quality}")
    elif encoding == "phred64":
        if not 0 <= quality <= 62:
            raise ValueError(f"Phred+64 quality must be 0-62, got {quality}")
    else:
        raise ValueError(f"Unknown encoding: {encoding}")


def validate_sequence(sequence: str, allowed_bases: str = "ACGTN") -> None:
    """
    Validate DNA/RNA sequence.

    Args:
        sequence: Sequence string to validate
        allowed_bases: String of allowed base characters

    Raises:
        ValueError: If sequence contains invalid characters
    """
    sequence = sequence.upper()
    allowed_set = set(allowed_bases.upper())

    invalid_bases = set(sequence) - allowed_set
    if invalid_bases:
        raise ValueError(
            f"Sequence contains invalid bases: {invalid_bases}. "
            f"Allowed: {allowed_bases}"
        )


def validate_adapter_sequence(adapter: str) -> None:
    """
    Validate adapter sequence.

    Args:
        adapter: Adapter sequence string

    Raises:
        ValueError: If adapter is invalid
    """
    if not adapter:
        raise ValueError("Adapter sequence cannot be empty")

    if len(adapter) < 6:
        raise ValueError("Adapter sequence too short (minimum 6 bases)")

    validate_sequence(adapter, allowed_bases="ACGT")


def validate_threads(threads: int) -> None:
    """
    Validate number of threads.

    Args:
        threads: Number of threads

    Raises:
        ValueError: If threads is invalid
    """
    if threads < 1:
        raise ValueError("Number of threads must be at least 1")

    if threads > 128:
        raise ValueError("Number of threads too high (maximum 128)")


def validate_output_directory(output_dir: Path | str, create: bool = True) -> Path:
    """
    Validate and optionally create output directory.

    Args:
        output_dir: Output directory path
        create: Whether to create directory if it doesn't exist

    Returns:
        Validated Path object

    Raises:
        ValueError: If path exists but is not a directory
    """
    output_dir = Path(output_dir)

    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"Path exists but is not a directory: {output_dir}")

    if create and not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir
