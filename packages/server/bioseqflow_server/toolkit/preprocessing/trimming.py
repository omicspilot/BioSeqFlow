from __future__ import annotations

"""Adapter trimming and quality trimming modules."""

import subprocess
from pathlib import Path
from typing import Any

from bioseqflow_server.toolkit.core.base import PreprocessingModule
from bioseqflow_server.toolkit.utils.alignment import find_adapter_fuzzy
from bioseqflow_server.toolkit.utils.io import FastqRecord, read_fastq, write_fastq
from bioseqflow_server.toolkit.utils.validators import (
    validate_adapter_sequence,
    validate_file_path,
    validate_quality_score,
)


class AdapterTrimmer(PreprocessingModule):
    """Trim adapter sequences from reads."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize adapter trimmer.

        Args:
            config: Configuration object
        """
        super().__init__(config)
        self.min_overlap = 3
        self.error_rate = 0.1

    def validate_inputs(
        self, input_file: Path | str, _output_file: Path | str, adapter: str
    ) -> None:
        """
        Validate inputs for trimming.

        Args:
            input_file: Input FASTQ file
            _output_file: Output FASTQ file (unused in base validation)
            adapter: Adapter sequence

        Raises:
            ValueError: If inputs are invalid
        """
        validate_file_path(
            input_file,
            must_exist=True,
            extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        )
        validate_adapter_sequence(adapter)

    def process(self, input_file: Path, output_file: Path) -> dict[str, int | float]:
        """
        Process file (alias for trim method).

        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file

        Returns:
            Trimming statistics
        """
        if self.config and hasattr(self.config, "adapter_sequence"):
            adapter = self.config.adapter_sequence
        else:
            raise ValueError("Adapter sequence not provided in config")

        return self.trim(input_file, output_file, adapter)

    def trim(
        self,
        input_file: Path | str,
        output_file: Path | str,
        adapter: str,
        use_cutadapt: bool = True,
        fuzzy_match: bool = False,
        max_error_rate: float = 0.15,
    ) -> dict[str, int | float]:
        """
        Trim adapter sequences from reads.

        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            adapter: Adapter sequence to trim
            use_cutadapt: Use cutadapt if available, else use built-in trimmer
            fuzzy_match: Use Smith-Waterman fuzzy matching (only for built-in)
            max_error_rate: Maximum error rate for fuzzy matching (0.0-1.0)

        Returns:
            Dictionary with trimming statistics

        Raises:
            RuntimeError: If trimming fails

        Notes:
            - Fuzzy matching uses Smith-Waterman local alignment
            - Allows mismatches and partial adapter matches
            - Recommended for degraded samples or high error rates
            - cutadapt always uses fuzzy matching internally
        """
        input_file = Path(input_file)
        output_file = Path(output_file)

        self.validate_inputs(input_file, output_file, adapter)

        if use_cutadapt:
            try:
                return self._trim_with_cutadapt(input_file, output_file, adapter)
            except (FileNotFoundError, RuntimeError):
                # Fall back to built-in trimmer
                return self._trim_builtin(
                    input_file, output_file, adapter, fuzzy_match, max_error_rate
                )
        else:
            return self._trim_builtin(input_file, output_file, adapter, fuzzy_match, max_error_rate)

    def _trim_with_cutadapt(
        self, input_file: Path, output_file: Path, adapter: str
    ) -> dict[str, int | float]:
        """
        Trim using cutadapt tool.

        Args:
            input_file: Input file
            output_file: Output file
            adapter: Adapter sequence

        Returns:
            Trimming statistics
        """
        output_file.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "cutadapt",
            "-a",
            adapter,
            "-o",
            str(output_file),
            str(input_file),
            "--overlap",
            str(self.min_overlap),
            "--error-rate",
            str(self.error_rate),
        ]

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Parse cutadapt output
            self.stats = self._parse_cutadapt_output(result.stdout)
            return self.stats

        except FileNotFoundError as e:
            raise RuntimeError("cutadapt not found. Install with: pip install cutadapt") from e

    def _trim_builtin(
        self,
        input_file: Path,
        output_file: Path,
        adapter: str,
        fuzzy_match: bool = False,
        max_error_rate: float = 0.15,
    ) -> dict[str, int | float]:
        """
        Built-in adapter trimming with optional fuzzy matching.

        Args:
            input_file: Input file
            output_file: Output file
            adapter: Adapter sequence
            fuzzy_match: Use Smith-Waterman fuzzy matching
            max_error_rate: Maximum error rate for fuzzy matching

        Returns:
            Trimming statistics
        """
        output_file.parent.mkdir(parents=True, exist_ok=True)

        total_reads = 0
        trimmed_reads = 0
        total_bp_removed = 0

        def trim_record(record: FastqRecord) -> FastqRecord:
            nonlocal total_reads, trimmed_reads, total_bp_removed
            total_reads += 1

            adapter_pos = -1

            if fuzzy_match:
                # Use Smith-Waterman fuzzy matching
                result = find_adapter_fuzzy(
                    record.sequence,
                    adapter,
                    min_overlap=self.min_overlap,
                    max_error_rate=max_error_rate,
                )
                if result is not None:
                    adapter_pos = result[0]  # Start position
            else:
                # Use exact string matching
                adapter_pos = record.sequence.find(adapter)

            if adapter_pos != -1:
                trimmed_reads += 1
                bp_removed = len(record.sequence) - adapter_pos
                total_bp_removed += bp_removed

                # Trim sequence and quality
                return FastqRecord(
                    record.header,
                    record.sequence[:adapter_pos],
                    record.plus,
                    record.quality[:adapter_pos],
                )

            return record

        # Process records
        trimmed_records = (trim_record(record) for record in read_fastq(input_file))
        write_fastq(trimmed_records, output_file, compress=str(output_file).endswith(".gz"))

        self.stats = {
            "total_reads": total_reads,
            "trimmed_reads": trimmed_reads,
            "percent_trimmed": (trimmed_reads / total_reads * 100) if total_reads > 0 else 0,
            "total_bp_removed": total_bp_removed,
        }

        return self.stats

    def _parse_cutadapt_output(self, stdout: str) -> dict[str, int | float]:
        """
        Parse cutadapt stdout for statistics.

        Args:
            stdout: Cutadapt stdout string

        Returns:
            Statistics dictionary
        """
        stats: dict[str, int | float] = {}

        for line in stdout.split("\n"):
            if "Total reads processed:" in line:
                stats["total_reads"] = int(line.split()[-1].replace(",", ""))
            elif "Reads with adapters:" in line:
                parts = line.split()
                stats["trimmed_reads"] = int(parts[-2].replace(",", ""))
            elif "Total basepairs processed:" in line:
                stats["total_bp"] = int(line.split()[-2].replace(",", ""))
            elif "Total written (filtered):" in line:
                stats["written_bp"] = int(line.split()[-3].replace(",", ""))

        # Calculate percentage
        if "total_reads" in stats and "trimmed_reads" in stats:
            stats["percent_trimmed"] = stats["trimmed_reads"] / stats["total_reads"] * 100

        return stats


class QualityTrimmer(PreprocessingModule):
    """Trim low-quality bases from read ends."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize quality trimmer.

        Args:
            config: Configuration object
        """
        super().__init__(config)

    def validate_inputs(
        self, input_file: Path | str, _output_file: Path | str, min_quality: int
    ) -> None:
        """
        Validate inputs.

        Args:
            input_file: Input file
            _output_file: Output file (unused in base validation)
            min_quality: Minimum quality threshold

        Raises:
            ValueError: If inputs are invalid
        """
        validate_file_path(
            input_file,
            must_exist=True,
            extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        )
        validate_quality_score(min_quality)

    def process(self, input_file: Path, output_file: Path) -> dict[str, int | float]:
        """
        Process file with quality trimming.

        Args:
            input_file: Input file
            output_file: Output file

        Returns:
            Processing statistics
        """
        min_quality = (
            self.config.min_quality if self.config and hasattr(self.config, "min_quality") else 20
        )
        return self.trim(input_file, output_file, min_quality)

    def trim(
        self,
        input_file: Path | str,
        output_file: Path | str,
        min_quality: int = 20,
        trim_5prime: bool = True,
        trim_3prime: bool = True,
    ) -> dict[str, int | float]:
        """
        Trim low-quality bases from read ends.

        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            min_quality: Minimum Phred quality score
            trim_5prime: Trim from 5' end
            trim_3prime: Trim from 3' end

        Returns:
            Trimming statistics
        """
        input_file = Path(input_file)
        output_file = Path(output_file)

        self.validate_inputs(input_file, output_file, min_quality)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        total_reads = 0
        trimmed_reads = 0
        total_bp_removed = 0

        def trim_record(record: FastqRecord) -> FastqRecord | None:
            nonlocal total_reads, trimmed_reads, total_bp_removed
            total_reads += 1

            start = 0
            end = len(record.quality)
            original_length = end

            # Trim from 5' end
            if trim_5prime:
                for i, q in enumerate(record.quality):
                    if ord(q) - 33 >= min_quality:
                        start = i
                        break
                else:
                    # No base met quality threshold from 5' end
                    start = original_length

            # Trim from 3' end
            if trim_3prime:
                for i in range(len(record.quality) - 1, -1, -1):
                    if ord(record.quality[i]) - 33 >= min_quality:
                        end = i + 1
                        break
                else:
                    # No base met quality threshold from 3' end
                    end = 0

            # CRITICAL FIX: Check if read is entirely low quality
            if start >= end:
                # Entire read is low quality - filter it out
                trimmed_reads += 1
                total_bp_removed += original_length
                return None

            # Check if trimming occurred
            if start > 0 or end < original_length:
                trimmed_reads += 1
                total_bp_removed += start + (original_length - end)

            return FastqRecord(
                record.header,
                record.sequence[start:end],
                record.plus,
                record.quality[start:end],
            )

        # Process records and filter out None (all-low-quality reads)
        trimmed_records = (
            trimmed
            for record in read_fastq(input_file)
            if (trimmed := trim_record(record)) is not None
        )
        written = write_fastq(
            trimmed_records, output_file, compress=str(output_file).endswith(".gz")
        )

        self.stats = {
            "total_reads": total_reads,
            "trimmed_reads": trimmed_reads,
            "reads_written": written,
            "reads_discarded": total_reads - written,
            "percent_trimmed": (trimmed_reads / total_reads * 100) if total_reads > 0 else 0,
            "total_bp_removed": total_bp_removed,
        }

        return self.stats
