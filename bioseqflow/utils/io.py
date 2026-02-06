from __future__ import annotations

"""I/O utilities for reading and writing sequence files."""

import gzip
from pathlib import Path
from typing import Generator, TextIO

import pandas as pd


class FastqRecord:
    """Represents a single FASTQ record."""

    def __init__(
        self, header: str, sequence: str, plus: str, quality: str
    ) -> None:
        """
        Initialize FASTQ record.

        Args:
            header: Sequence identifier line (starts with @)
            sequence: DNA/RNA sequence
            plus: Separator line (starts with +)
            quality: Quality scores (Phred+33 encoded)
        """
        self.header = header.strip()
        self.sequence = sequence.strip()
        self.plus = plus.strip()
        self.quality = quality.strip()

    @property
    def name(self) -> str:
        """Get sequence name without @ prefix."""
        return self.header[1:] if self.header.startswith("@") else self.header

    @property
    def length(self) -> int:
        """Get sequence length."""
        return len(self.sequence)

    def mean_quality(self) -> float:
        """
        Calculate mean quality score.

        Returns:
            Mean Phred quality score
        """
        if not self.quality:
            return 0.0
        return sum(ord(q) - 33 for q in self.quality) / len(self.quality)

    def __str__(self) -> str:
        """String representation (4-line FASTQ format)."""
        return f"{self.header}\n{self.sequence}\n{self.plus}\n{self.quality}\n"


def open_file(file_path: Path | str, mode: str = "r") -> TextIO:
    """
    Open a file, handling gzip compression automatically.

    Args:
        file_path: Path to file
        mode: File mode ('r' or 'w')

    Returns:
        File handle

    Raises:
        FileNotFoundError: If file doesn't exist (read mode)
    """
    file_path = Path(file_path)

    if "r" in mode and not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    if file_path.suffix == ".gz":
        return gzip.open(file_path, mode + "t", encoding="utf-8")
    return open(file_path, mode, encoding="utf-8")


def read_fastq(file_path: Path | str) -> Generator[FastqRecord, None, None]:
    """
    Read FASTQ file and yield records.

    Args:
        file_path: Path to FASTQ file (.fastq or .fastq.gz)

    Yields:
        FastqRecord objects

    Raises:
        ValueError: If FASTQ format is invalid
    """
    with open_file(file_path, "r") as f:
        line_num = 0
        while True:
            header = f.readline()
            line_num += 1
            if not header:
                break

            sequence = f.readline()
            line_num += 1
            plus = f.readline()
            line_num += 1
            quality = f.readline()
            line_num += 1

            if not (header and sequence and plus and quality):
                raise ValueError(
                    f"Incomplete FASTQ record at line {line_num - 3}"
                )

            if not header.startswith("@"):
                raise ValueError(
                    f"Invalid FASTQ header at line {line_num - 3}: {header.strip()}"
                )

            # CRITICAL FIX: Validate sequence and quality lengths match
            seq_len = len(sequence.strip())
            qual_len = len(quality.strip())
            if seq_len != qual_len:
                raise ValueError(
                    f"Sequence/quality length mismatch at line {line_num - 3} "
                    f"(header: {header.strip()}): "
                    f"sequence={seq_len}bp, quality={qual_len}bp"
                )

            yield FastqRecord(header, sequence, plus, quality)


def write_fastq(
    records: list[FastqRecord] | Generator[FastqRecord, None, None],
    file_path: Path | str,
    compress: bool = False,
) -> int:
    """
    Write FASTQ records to file.

    Args:
        records: Iterable of FastqRecord objects
        file_path: Output file path
        compress: Whether to gzip compress output

    Returns:
        Number of records written
    """
    file_path = Path(file_path)
    if compress and not str(file_path).endswith(".gz"):
        file_path = Path(str(file_path) + ".gz")

    file_path.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with open_file(file_path, "w") as f:
        for record in records:
            f.write(str(record))
            count += 1

    return count


def parse_sample_sheet(
    csv_path: Path | str,
) -> list[dict[str, str]]:
    """
    Parse sample sheet CSV file.

    Expected columns:
    - sample_id: Unique sample identifier
    - fastq_r1: Path to R1 FASTQ file
    - fastq_r2: Path to R2 FASTQ file (optional)
    - adapter: Adapter sequence (optional)
    - min_quality: Minimum quality threshold (optional)

    Args:
        csv_path: Path to CSV sample sheet

    Returns:
        List of sample dictionaries

    Raises:
        FileNotFoundError: If CSV file doesn't exist
        ValueError: If required columns are missing
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {csv_path}")

    df = pd.read_csv(csv_path)

    required_cols = ["sample_id", "fastq_r1"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    return df.to_dict("records")


def count_reads(file_path: Path | str) -> int:
    """
    Count number of reads in FASTQ file.

    Args:
        file_path: Path to FASTQ file

    Returns:
        Number of reads
    """
    count = 0
    for _ in read_fastq(file_path):
        count += 1
    return count


def get_read_length_distribution(
    file_path: Path | str, sample_size: int | None = None
) -> dict[int, int]:
    """
    Get read length distribution.

    Args:
        file_path: Path to FASTQ file
        sample_size: Number of reads to sample (None = all)

    Returns:
        Dictionary mapping length to count
    """
    length_dist: dict[int, int] = {}
    count = 0

    for record in read_fastq(file_path):
        length = record.length
        length_dist[length] = length_dist.get(length, 0) + 1
        count += 1

        if sample_size and count >= sample_size:
            break

    return length_dist
