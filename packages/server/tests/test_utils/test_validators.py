"""Tests for validation utilities."""

from pathlib import Path

import pytest

from bioseqflow_server.toolkit.utils.validators import (
    validate_adapter_sequence,
    validate_file_path,
    validate_output_directory,
    validate_quality_score,
    validate_sequence,
    validate_threads,
)


class TestFilePathValidation:
    """Tests for file path validation."""

    def test_validate_existing_file(self, sample_fastq_file: Path) -> None:
        """Test validating existing file."""
        result = validate_file_path(sample_fastq_file, must_exist=True)

        assert result == sample_fastq_file
        assert result.exists()

    def test_validate_nonexistent_file_required(self, temp_dir: Path) -> None:
        """Test validating non-existent file when required."""
        with pytest.raises(FileNotFoundError):
            validate_file_path(temp_dir / "nonexistent.fastq", must_exist=True)

    def test_validate_nonexistent_file_optional(self, temp_dir: Path) -> None:
        """Test validating non-existent file when optional."""
        result = validate_file_path(temp_dir / "new_file.fastq", must_exist=False)

        assert isinstance(result, Path)

    def test_validate_file_extension(self, sample_fastq_file: Path) -> None:
        """Test validating file extension."""
        result = validate_file_path(
            sample_fastq_file,
            must_exist=True,
            extensions=[".fastq", ".fq", ".fastq.gz"],
        )

        assert result == sample_fastq_file

    def test_validate_wrong_extension(self, temp_dir: Path) -> None:
        """Test validating wrong file extension."""
        wrong_file = temp_dir / "test.txt"
        wrong_file.write_text("content")

        with pytest.raises(ValueError, match="Invalid file extension"):
            validate_file_path(wrong_file, must_exist=True, extensions=[".fastq"])


class TestQualityScoreValidation:
    """Tests for quality score validation."""

    def test_validate_phred33_valid(self) -> None:
        """Test validating valid Phred+33 scores."""
        validate_quality_score(0, encoding="phred33")
        validate_quality_score(30, encoding="phred33")
        validate_quality_score(93, encoding="phred33")

    def test_validate_phred33_invalid(self) -> None:
        """Test validating invalid Phred+33 scores."""
        with pytest.raises(ValueError):
            validate_quality_score(-1, encoding="phred33")

        with pytest.raises(ValueError):
            validate_quality_score(94, encoding="phred33")

    def test_validate_unknown_encoding(self) -> None:
        """Test validating with unknown encoding."""
        with pytest.raises(ValueError, match="Unknown encoding"):
            validate_quality_score(30, encoding="unknown")


class TestSequenceValidation:
    """Tests for sequence validation."""

    def test_validate_dna_sequence(self) -> None:
        """Test validating valid DNA sequence."""
        validate_sequence("ATCGATCG")
        validate_sequence("atcgatcg")  # Should handle lowercase
        validate_sequence("ATCGNNNATCG")  # N is allowed

    def test_validate_invalid_sequence(self) -> None:
        """Test validating invalid sequence."""
        with pytest.raises(ValueError, match="invalid bases"):
            validate_sequence("ATCGXYZ")


class TestAdapterValidation:
    """Tests for adapter sequence validation."""

    def test_validate_valid_adapter(self) -> None:
        """Test validating valid adapter."""
        validate_adapter_sequence("AGATCGGAAGAG")
        validate_adapter_sequence("ATCGATCG")

    def test_validate_empty_adapter(self) -> None:
        """Test validating empty adapter."""
        with pytest.raises(ValueError, match="cannot be empty"):
            validate_adapter_sequence("")

    def test_validate_short_adapter(self) -> None:
        """Test validating too short adapter."""
        with pytest.raises(ValueError, match="too short"):
            validate_adapter_sequence("ATCG")

    def test_validate_adapter_with_invalid_bases(self) -> None:
        """Test validating adapter with invalid bases."""
        with pytest.raises(ValueError):
            validate_adapter_sequence("ATCGNN")  # N not allowed in adapters


class TestThreadsValidation:
    """Tests for threads validation."""

    def test_validate_valid_threads(self) -> None:
        """Test validating valid thread counts."""
        validate_threads(1)
        validate_threads(4)
        validate_threads(64)

    def test_validate_zero_threads(self) -> None:
        """Test validating zero threads."""
        with pytest.raises(ValueError, match="at least 1"):
            validate_threads(0)

    def test_validate_negative_threads(self) -> None:
        """Test validating negative threads."""
        with pytest.raises(ValueError):
            validate_threads(-1)

    def test_validate_too_many_threads(self) -> None:
        """Test validating excessive threads."""
        with pytest.raises(ValueError, match="too high"):
            validate_threads(200)


class TestOutputDirectoryValidation:
    """Tests for output directory validation."""

    def test_validate_existing_directory(self, temp_dir: Path) -> None:
        """Test validating existing directory."""
        result = validate_output_directory(temp_dir, create=False)

        assert result == temp_dir
        assert result.is_dir()

    def test_validate_create_directory(self, temp_dir: Path) -> None:
        """Test creating new directory."""
        new_dir = temp_dir / "new_output"
        result = validate_output_directory(new_dir, create=True)

        assert result.exists()
        assert result.is_dir()

    def test_validate_file_as_directory(self, temp_dir: Path) -> None:
        """Test validating file path as directory."""
        file_path = temp_dir / "file.txt"
        file_path.write_text("content")

        with pytest.raises(ValueError, match="not a directory"):
            validate_output_directory(file_path, create=False)
