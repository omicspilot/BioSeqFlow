from __future__ import annotations

"""Duplicate removal modules."""

from collections import defaultdict
from pathlib import Path
from typing import Any

from bioseqflow.core.base import PreprocessingModule
from bioseqflow.utils.io import FastqRecord, read_fastq, write_fastq
from bioseqflow.utils.validators import validate_file_path


class DuplicateRemover(PreprocessingModule):
    """Remove duplicate reads."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize duplicate remover.

        Args:
            config: Configuration object
        """
        super().__init__(config)

    def validate_inputs(
        self, input_file: Path | str, output_file: Path | str
    ) -> None:
        """
        Validate inputs.

        Args:
            input_file: Input file
            output_file: Output file

        Raises:
            ValueError: If inputs are invalid
        """
        validate_file_path(
            input_file,
            must_exist=True,
            extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        )

    def process(self, input_file: Path, output_file: Path) -> dict[str, int | float]:
        """
        Process file with duplicate removal.

        Args:
            input_file: Input file
            output_file: Output file

        Returns:
            Deduplication statistics
        """
        return self.remove_duplicates(input_file, output_file)

    def remove_duplicates(
        self,
        input_file: Path | str,
        output_file: Path | str,
        method: str = "exact",
        keep_best: bool = True,
    ) -> dict[str, int | float]:
        """
        Remove duplicate reads.

        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            method: Deduplication method ('exact' or 'prefix')
            keep_best: Keep read with highest quality (if duplicates found)

        Returns:
            Deduplication statistics

        Raises:
            ValueError: If method is invalid
        """
        input_file = Path(input_file)
        output_file = Path(output_file)

        self.validate_inputs(input_file, output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        if method == "exact":
            return self._remove_exact_duplicates(input_file, output_file, keep_best)
        elif method == "prefix":
            return self._remove_prefix_duplicates(input_file, output_file, keep_best)
        else:
            raise ValueError(f"Unknown deduplication method: {method}")

    def _remove_exact_duplicates(
        self, input_file: Path, output_file: Path, keep_best: bool
    ) -> dict[str, int | float]:
        """
        Remove exact sequence duplicates.

        Args:
            input_file: Input file
            output_file: Output file
            keep_best: Keep highest quality duplicate

        Returns:
            Statistics
        """
        total_reads = 0
        unique_reads = 0
        duplicates = 0

        seen_sequences: dict[str, FastqRecord] = {}

        for record in read_fastq(input_file):
            total_reads += 1

            if record.sequence not in seen_sequences:
                seen_sequences[record.sequence] = record
                unique_reads += 1
            else:
                duplicates += 1

                # Keep record with better quality if requested
                if keep_best:
                    if record.mean_quality() > seen_sequences[record.sequence].mean_quality():
                        seen_sequences[record.sequence] = record

        # Write unique reads
        write_fastq(
            seen_sequences.values(),
            output_file,
            compress=str(output_file).endswith(".gz")
        )

        self.stats = {
            "total_reads": total_reads,
            "unique_reads": unique_reads,
            "duplicates": duplicates,
            "duplication_rate": (duplicates / total_reads * 100) if total_reads > 0 else 0,
        }

        return self.stats

    def _remove_prefix_duplicates(
        self, input_file: Path, output_file: Path, keep_best: bool, prefix_length: int = 30
    ) -> dict[str, int | float]:
        """
        Remove duplicates based on sequence prefix.

        Useful for identifying PCR duplicates where only the start is identical.

        Args:
            input_file: Input file
            output_file: Output file
            keep_best: Keep highest quality duplicate
            prefix_length: Length of prefix to compare

        Returns:
            Statistics
        """
        total_reads = 0
        unique_reads = 0
        duplicates = 0

        seen_prefixes: dict[str, FastqRecord] = {}

        for record in read_fastq(input_file):
            total_reads += 1

            # Get prefix
            prefix = record.sequence[:prefix_length]

            if prefix not in seen_prefixes:
                seen_prefixes[prefix] = record
                unique_reads += 1
            else:
                duplicates += 1

                # Keep record with better quality if requested
                if keep_best:
                    if record.mean_quality() > seen_prefixes[prefix].mean_quality():
                        seen_prefixes[prefix] = record

        # Write unique reads
        write_fastq(
            seen_prefixes.values(),
            output_file,
            compress=str(output_file).endswith(".gz")
        )

        self.stats = {
            "total_reads": total_reads,
            "unique_reads": unique_reads,
            "duplicates": duplicates,
            "duplication_rate": (duplicates / total_reads * 100) if total_reads > 0 else 0,
            "method": "prefix",
            "prefix_length": prefix_length,
        }

        return self.stats


class UMIDeduplicator(PreprocessingModule):
    """Remove duplicates using Unique Molecular Identifiers (UMIs)."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize UMI deduplicator.

        Args:
            config: Configuration object
        """
        super().__init__(config)

    def validate_inputs(
        self, input_file: Path | str, output_file: Path | str
    ) -> None:
        """
        Validate inputs.

        Args:
            input_file: Input file
            output_file: Output file

        Raises:
            ValueError: If inputs are invalid
        """
        validate_file_path(
            input_file,
            must_exist=True,
            extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        )

    def process(self, input_file: Path, output_file: Path) -> dict[str, int | float]:
        """
        Process file with UMI-based deduplication.

        Args:
            input_file: Input file
            output_file: Output file

        Returns:
            Deduplication statistics
        """
        return self.deduplicate(input_file, output_file)

    def deduplicate(
        self,
        input_file: Path | str,
        output_file: Path | str,
        umi_length: int = 8,
        umi_location: str = "start",
    ) -> dict[str, int | float]:
        """
        Deduplicate reads using UMIs.

        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            umi_length: Length of UMI sequence
            umi_location: UMI location ('start' or 'end')

        Returns:
            Deduplication statistics
        """
        input_file = Path(input_file)
        output_file = Path(output_file)

        self.validate_inputs(input_file, output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        total_reads = 0
        unique_umis = 0
        duplicates = 0

        # Group by UMI
        umi_groups: dict[str, list[FastqRecord]] = defaultdict(list)

        for record in read_fastq(input_file):
            total_reads += 1

            # Extract UMI
            if umi_location == "start":
                umi = record.sequence[:umi_length]
                trimmed_seq = record.sequence[umi_length:]
                trimmed_qual = record.quality[umi_length:]
            else:  # end
                umi = record.sequence[-umi_length:]
                trimmed_seq = record.sequence[:-umi_length]
                trimmed_qual = record.quality[:-umi_length]

            # Create trimmed record (without UMI)
            trimmed_record = FastqRecord(
                record.header, trimmed_seq, record.plus, trimmed_qual
            )

            umi_groups[umi].append(trimmed_record)

        # Keep best quality read from each UMI group
        unique_records = []
        for umi, records in umi_groups.items():
            unique_umis += 1
            duplicates += len(records) - 1

            # Keep read with highest mean quality
            best_record = max(records, key=lambda r: r.mean_quality())
            unique_records.append(best_record)

        # Write deduplicated reads
        write_fastq(unique_records, output_file, compress=str(output_file).endswith(".gz"))

        self.stats = {
            "total_reads": total_reads,
            "unique_umis": unique_umis,
            "duplicates": duplicates,
            "duplication_rate": (duplicates / total_reads * 100) if total_reads > 0 else 0,
            "umi_length": umi_length,
        }

        return self.stats
