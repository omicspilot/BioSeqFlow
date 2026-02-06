"""Tests for filtering modules."""

from pathlib import Path

import pytest

from bioseqflow.preprocessing.filtering import LengthFilter, QualityFilter


class TestQualityFilter:
    """Tests for QualityFilter class."""

    @pytest.fixture
    def filter_module(self) -> QualityFilter:
        """Create QualityFilter instance."""
        return QualityFilter()

    def test_filter_by_quality(
        self, filter_module: QualityFilter, sample_fastq_file: Path, temp_dir: Path
    ) -> None:
        """Test filtering by quality score."""
        output_file = temp_dir / "filtered.fastq"

        stats = filter_module.filter(
            sample_fastq_file,
            output_file,
            min_quality=10,  # Low threshold to pass some reads
            min_length=10,
        )

        assert "total_reads" in stats
        assert "passed_reads" in stats
        assert "percent_passed" in stats
        assert output_file.exists()

    def test_filter_by_length(
        self, filter_module: QualityFilter, sample_fastq_file: Path, temp_dir: Path
    ) -> None:
        """Test filtering by read length."""
        output_file = temp_dir / "filtered.fastq"

        stats = filter_module.filter(
            sample_fastq_file, output_file, min_quality=0, min_length=100
        )

        assert stats["total_reads"] >= 0
        assert stats["passed_reads"] <= stats["total_reads"]

    def test_filter_statistics(
        self, filter_module: QualityFilter, sample_fastq_file: Path, temp_dir: Path
    ) -> None:
        """Test filtering statistics."""
        output_file = temp_dir / "filtered.fastq"

        stats = filter_module.filter(
            sample_fastq_file, output_file, min_quality=0, min_length=10
        )

        assert "failed_quality" in stats
        assert "failed_length" in stats
        assert "failed_n_content" in stats
        assert (
            stats["passed_reads"]
            + stats["failed_quality"]
            + stats["failed_length"]
            + stats["failed_n_content"]
            <= stats["total_reads"]
        )


class TestLengthFilter:
    """Tests for LengthFilter class."""

    @pytest.fixture
    def filter_module(self) -> LengthFilter:
        """Create LengthFilter instance."""
        return LengthFilter()

    def test_filter_min_length(
        self, filter_module: LengthFilter, sample_fastq_file: Path, temp_dir: Path
    ) -> None:
        """Test filtering by minimum length."""
        output_file = temp_dir / "length_filtered.fastq"

        stats = filter_module.filter(sample_fastq_file, output_file, min_length=50)

        assert output_file.exists()
        assert stats["total_reads"] >= 0

    def test_filter_max_length(
        self, filter_module: LengthFilter, sample_fastq_file: Path, temp_dir: Path
    ) -> None:
        """Test filtering by maximum length."""
        output_file = temp_dir / "length_filtered.fastq"

        stats = filter_module.filter(
            sample_fastq_file, output_file, min_length=10, max_length=100
        )

        assert "too_long" in stats
        assert "too_short" in stats

    def test_filter_length_statistics(
        self, filter_module: LengthFilter, sample_fastq_file: Path, temp_dir: Path
    ) -> None:
        """Test length filter statistics."""
        output_file = temp_dir / "length_filtered.fastq"

        stats = filter_module.filter(
            sample_fastq_file, output_file, min_length=30, max_length=200
        )

        assert stats["total_reads"] == (
            stats["passed_reads"] + stats["too_short"] + stats["too_long"]
        )
