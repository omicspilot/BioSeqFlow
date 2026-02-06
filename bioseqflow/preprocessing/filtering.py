from __future__ import annotations

"""Quality and length filtering modules."""

from pathlib import Path
from typing import Any

from bioseqflow.core.base import PreprocessingModule
from bioseqflow.utils.io import FastqRecord, read_fastq, write_fastq
from bioseqflow.utils.validators import validate_file_path, validate_quality_score


class QualityFilter(PreprocessingModule):
    """Filter reads based on quality metrics."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize quality filter.

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
        Process file with quality filtering.

        Args:
            input_file: Input file
            output_file: Output file

        Returns:
            Filtering statistics
        """
        min_quality = (
            self.config.min_quality if self.config and hasattr(self.config, "min_quality") else 20
        )
        min_length = (
            self.config.min_length if self.config and hasattr(self.config, "min_length") else 50
        )

        return self.filter(input_file, output_file, min_quality=min_quality, min_length=min_length)

    def filter(
        self,
        input_file: Path | str,
        output_file: Path | str,
        min_quality: int = 20,
        min_length: int = 50,
        max_n_content: float = 5.0,
    ) -> dict[str, int | float]:
        """
        Filter reads based on quality and length.

        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            min_quality: Minimum mean quality score
            min_length: Minimum read length
            max_n_content: Maximum N content percentage

        Returns:
            Filtering statistics
        """
        input_file = Path(input_file)
        output_file = Path(output_file)

        self.validate_inputs(input_file, output_file)
        validate_quality_score(min_quality)

        output_file.parent.mkdir(parents=True, exist_ok=True)

        total_reads = 0
        passed_reads = 0
        failed_quality = 0
        failed_length = 0
        failed_n_content = 0

        def passes_filters(record: FastqRecord) -> bool:
            nonlocal total_reads, passed_reads, failed_quality, failed_length, failed_n_content
            total_reads += 1

            # OPTIMIZATION: Check in order of computational cost and rejection rate
            # 1. Mean quality (O(n) but fast, high rejection rate)
            if record.mean_quality() < min_quality:
                failed_quality += 1
                return False

            # 2. N content (O(n) string count, medium rejection rate)
            n_content = (record.sequence.count("N") / record.length * 100) if record.length > 0 else 0
            if n_content > max_n_content:
                failed_n_content += 1
                return False

            # 3. Length (O(1), but less discriminatory)
            if record.length < min_length:
                failed_length += 1
                return False

            passed_reads += 1
            return True

        # Filter records
        filtered_records = (record for record in read_fastq(input_file) if passes_filters(record))
        write_fastq(filtered_records, output_file, compress=str(output_file).endswith(".gz"))

        self.stats = {
            "total_reads": total_reads,
            "passed_reads": passed_reads,
            "failed_quality": failed_quality,
            "failed_length": failed_length,
            "failed_n_content": failed_n_content,
            "percent_passed": (passed_reads / total_reads * 100) if total_reads > 0 else 0,
        }

        return self.stats


class LengthFilter(PreprocessingModule):
    """Filter reads based on length."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize length filter.

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
        Process file with length filtering.

        Args:
            input_file: Input file
            output_file: Output file

        Returns:
            Filtering statistics
        """
        min_length = (
            self.config.min_length if self.config and hasattr(self.config, "min_length") else 50
        )
        max_length = (
            self.config.max_length if self.config and hasattr(self.config, "max_length") else None
        )

        return self.filter(input_file, output_file, min_length=min_length, max_length=max_length)

    def filter(
        self,
        input_file: Path | str,
        output_file: Path | str,
        min_length: int = 50,
        max_length: int | None = None,
    ) -> dict[str, int | float]:
        """
        Filter reads based on length.

        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            min_length: Minimum read length
            max_length: Maximum read length (None = no maximum)

        Returns:
            Filtering statistics
        """
        input_file = Path(input_file)
        output_file = Path(output_file)

        self.validate_inputs(input_file, output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        total_reads = 0
        passed_reads = 0
        too_short = 0
        too_long = 0

        def passes_length_filter(record: FastqRecord) -> bool:
            nonlocal total_reads, passed_reads, too_short, too_long
            total_reads += 1

            if record.length < min_length:
                too_short += 1
                return False

            if max_length is not None and record.length > max_length:
                too_long += 1
                return False

            passed_reads += 1
            return True

        # Filter records
        filtered_records = (
            record for record in read_fastq(input_file) if passes_length_filter(record)
        )
        write_fastq(filtered_records, output_file, compress=str(output_file).endswith(".gz"))

        self.stats = {
            "total_reads": total_reads,
            "passed_reads": passed_reads,
            "too_short": too_short,
            "too_long": too_long,
            "percent_passed": (passed_reads / total_reads * 100) if total_reads > 0 else 0,
        }

        return self.stats
