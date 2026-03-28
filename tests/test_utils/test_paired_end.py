from __future__ import annotations

"""Tests for paired-end sequencing utilities."""


from typing import TYPE_CHECKING

import pytest

from bioseqflow.utils.io import FastqRecord, write_fastq
from bioseqflow.utils.paired_end import (
    calculate_insert_size_distribution,
    check_read_orientation,
    extract_read_id,
    read_paired_fastq,
    validate_paired_files,
)

if TYPE_CHECKING:
    from pathlib import Path


class TestExtractReadID:
    """Test read ID extraction."""

    def test_simple_header(self):
        """Test simple read header."""
        header = "@read1"
        assert extract_read_id(header) == "read1"

    def test_illumina_old_format(self):
        """Test old Illumina format with /1 /2."""
        header = "@instrument:run:flowcell:lane:tile:x:y/1"
        expected = "instrument:run:flowcell:lane:tile:x:y"
        assert extract_read_id(header) == expected

        header2 = "@instrument:run:flowcell:lane:tile:x:y/2"
        assert extract_read_id(header2) == expected

    def test_illumina_new_format(self):
        """Test new Illumina format with space-separated info."""
        header = "@instrument:run:flowcell:lane:tile:x:y 1:N:0:ATCG"
        expected = "instrument:run:flowcell:lane:tile:x:y"
        assert extract_read_id(header) == expected

        header2 = "@instrument:run:flowcell:lane:tile:x:y 2:N:0:ATCG"
        assert extract_read_id(header2) == expected

    def test_sra_format(self):
        """Test SRA format."""
        header = "@SRR123456.1 1 length=150"
        assert extract_read_id(header) == "SRR123456.1"

        header2 = "@SRR123456.1 2 length=150"
        assert extract_read_id(header2) == "SRR123456.1"

    def test_with_description(self):
        """Test header with description."""
        header = "@read1 description here"
        assert extract_read_id(header) == "read1"


class TestValidatePairedFiles:
    """Test paired-end file validation."""

    def test_valid_paired_files(self, tmp_path: Path):
        """Test validation of properly paired files."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        # Create matching paired reads
        r1_records = [
            FastqRecord("@read1/1\n", "ATCG\n", "+\n", "IIII\n"),
            FastqRecord("@read2/1\n", "GCTA\n", "+\n", "IIII\n"),
            FastqRecord("@read3/1\n", "TTAA\n", "+\n", "IIII\n"),
        ]

        r2_records = [
            FastqRecord("@read1/2\n", "CGAT\n", "+\n", "IIII\n"),
            FastqRecord("@read2/2\n", "TAGC\n", "+\n", "IIII\n"),
            FastqRecord("@read3/2\n", "AATT\n", "+\n", "IIII\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        result = validate_paired_files(r1_file, r2_file)

        assert result.is_valid
        assert result.total_pairs == 3
        assert len(result.errors) == 0

    def test_mismatched_read_counts(self, tmp_path: Path):
        """Test validation fails with different read counts."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        # R1 has 3 reads, R2 has 2
        r1_records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
            FastqRecord("@read2\n", "GCTA\n", "+\n", "IIII\n"),
            FastqRecord("@read3\n", "TTAA\n", "+\n", "IIII\n"),
        ]

        r2_records = [
            FastqRecord("@read1\n", "CGAT\n", "+\n", "IIII\n"),
            FastqRecord("@read2\n", "TAGC\n", "+\n", "IIII\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        result = validate_paired_files(r1_file, r2_file)

        assert not result.is_valid
        assert "R2 file ended" in result.errors[0]

    def test_mismatched_read_ids(self, tmp_path: Path):
        """Test validation fails with mismatched read IDs."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        r1_records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
            FastqRecord("@read2\n", "GCTA\n", "+\n", "IIII\n"),
        ]

        r2_records = [
            FastqRecord("@read1\n", "CGAT\n", "+\n", "IIII\n"),
            FastqRecord("@read3\n", "TAGC\n", "+\n", "IIII\n"),  # Mismatched!
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        result = validate_paired_files(r1_file, r2_file)

        assert not result.is_valid
        assert "Read ID mismatch" in result.errors[0]

    def test_same_file_error(self, tmp_path: Path):
        """Test validation fails when R1 and R2 are same file."""
        r1_file = tmp_path / "test.fastq"

        records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
        ]
        write_fastq(records, r1_file)

        result = validate_paired_files(r1_file, r1_file)

        assert not result.is_valid
        assert "same file" in result.errors[0]

    def test_missing_file(self, tmp_path: Path):
        """Test validation fails with missing file."""
        r1_file = tmp_path / "exists.fastq"
        r2_file = tmp_path / "missing.fastq"

        records = [FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n")]
        write_fastq(records, r1_file)

        result = validate_paired_files(r1_file, r2_file)

        assert not result.is_valid
        assert "not found" in result.errors[0]

    def test_max_pairs_limit(self, tmp_path: Path):
        """Test validation with max pairs limit."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        # Create 10 pairs
        r1_records = [FastqRecord(f"@read{i}\n", "ATCG\n", "+\n", "IIII\n") for i in range(10)]
        r2_records = [FastqRecord(f"@read{i}\n", "CGAT\n", "+\n", "IIII\n") for i in range(10)]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        # Only check first 5 pairs
        result = validate_paired_files(r1_file, r2_file, max_pairs_to_check=5)

        assert result.is_valid
        assert result.total_pairs == 5

    def test_skip_order_check(self, tmp_path: Path):
        """Test validation without order checking."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        # Different read IDs but same count
        r1_records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
            FastqRecord("@read2\n", "GCTA\n", "+\n", "IIII\n"),
        ]

        r2_records = [
            FastqRecord("@read3\n", "CGAT\n", "+\n", "IIII\n"),
            FastqRecord("@read4\n", "TAGC\n", "+\n", "IIII\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        # Should pass without order check
        result = validate_paired_files(r1_file, r2_file, check_order=False)
        assert result.is_valid


class TestReadPairedFastq:
    """Test paired FASTQ reading."""

    def test_read_paired_basic(self, tmp_path: Path):
        """Test basic paired reading."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        r1_records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
            FastqRecord("@read2\n", "GCTA\n", "+\n", "IIII\n"),
        ]

        r2_records = [
            FastqRecord("@read1\n", "CGAT\n", "+\n", "IIII\n"),
            FastqRecord("@read2\n", "TAGC\n", "+\n", "IIII\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        pairs = list(read_paired_fastq(r1_file, r2_file))

        assert len(pairs) == 2
        assert pairs[0].r1.sequence.strip() == "ATCG"
        assert pairs[0].r2.sequence.strip() == "CGAT"
        assert pairs[1].r1.sequence.strip() == "GCTA"
        assert pairs[1].r2.sequence.strip() == "TAGC"

    def test_read_paired_mismatch_error(self, tmp_path: Path):
        """Test error on read ID mismatch."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        r1_records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
        ]

        r2_records = [
            FastqRecord("@read2\n", "CGAT\n", "+\n", "IIII\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        with pytest.raises(ValueError, match="Read ID mismatch"):
            list(read_paired_fastq(r1_file, r2_file, validate=True))

    def test_read_paired_unequal_length_error(self, tmp_path: Path):
        """Test error when files have different lengths."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        r1_records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
            FastqRecord("@read2\n", "GCTA\n", "+\n", "IIII\n"),
        ]

        r2_records = [
            FastqRecord("@read1\n", "CGAT\n", "+\n", "IIII\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        with pytest.raises(ValueError, match="R2 file ended early"):
            list(read_paired_fastq(r1_file, r2_file, validate=True))

    def test_read_paired_without_validation(self, tmp_path: Path):
        """Test reading without validation (faster)."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        r1_records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
        ]

        r2_records = [
            FastqRecord("@read2\n", "CGAT\n", "+\n", "IIII\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        # Should work without validation
        pairs = list(read_paired_fastq(r1_file, r2_file, validate=False))
        assert len(pairs) == 1


class TestInsertSizeDistribution:
    """Test insert size calculation."""

    def test_insert_size_basic(self, tmp_path: Path):
        """Test basic insert size calculation."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        # R1: 100bp, R2: 100bp -> estimated insert 200bp
        r1_records = [
            FastqRecord("@read1\n", "A" * 100 + "\n", "+\n", "I" * 100 + "\n"),
            FastqRecord("@read2\n", "A" * 100 + "\n", "+\n", "I" * 100 + "\n"),
        ]

        r2_records = [
            FastqRecord("@read1\n", "T" * 100 + "\n", "+\n", "I" * 100 + "\n"),
            FastqRecord("@read2\n", "T" * 100 + "\n", "+\n", "I" * 100 + "\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        stats = calculate_insert_size_distribution(r1_file, r2_file)

        assert stats["pairs_analyzed"] == 2
        assert stats["mean_insert_size"] == 200.0
        assert stats["min_insert_size"] == 200
        assert stats["max_insert_size"] == 200

    def test_insert_size_max_reads(self, tmp_path: Path):
        """Test insert size with max reads limit."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        # Create 10 pairs
        r1_records = [
            FastqRecord(f"@read{i}\n", "A" * 100 + "\n", "+\n", "I" * 100 + "\n") for i in range(10)
        ]
        r2_records = [
            FastqRecord(f"@read{i}\n", "T" * 100 + "\n", "+\n", "I" * 100 + "\n") for i in range(10)
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        # Only analyze 5 pairs
        stats = calculate_insert_size_distribution(r1_file, r2_file, max_reads=5)

        assert stats["pairs_analyzed"] == 5


class TestReadOrientation:
    """Test read orientation checking."""

    def test_consistent_lengths(self, tmp_path: Path):
        """Test reads with consistent lengths."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        # Both 100bp
        r1_records = [
            FastqRecord("@read1\n", "A" * 100 + "\n", "+\n", "I" * 100 + "\n"),
            FastqRecord("@read2\n", "A" * 100 + "\n", "+\n", "I" * 100 + "\n"),
        ]

        r2_records = [
            FastqRecord("@read1\n", "T" * 100 + "\n", "+\n", "I" * 100 + "\n"),
            FastqRecord("@read2\n", "T" * 100 + "\n", "+\n", "I" * 100 + "\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        stats = check_read_orientation(r1_file, r2_file)

        assert stats["pairs_checked"] == 2
        assert stats["mean_length_diff"] == 0.0
        assert stats["max_length_diff"] == 0
        assert stats["length_consistent"]

    def test_inconsistent_lengths(self, tmp_path: Path):
        """Test reads with inconsistent lengths."""
        r1_file = tmp_path / "test_R1.fastq"
        r2_file = tmp_path / "test_R2.fastq"

        # R1: 100bp, R2: 50bp (large difference)
        r1_records = [
            FastqRecord("@read1\n", "A" * 100 + "\n", "+\n", "I" * 100 + "\n"),
        ]

        r2_records = [
            FastqRecord("@read1\n", "T" * 50 + "\n", "+\n", "I" * 50 + "\n"),
        ]

        write_fastq(r1_records, r1_file)
        write_fastq(r2_records, r2_file)

        stats = check_read_orientation(r1_file, r2_file)

        assert stats["pairs_checked"] == 1
        assert stats["mean_length_diff"] == 50.0
        assert stats["max_length_diff"] == 50
        assert not stats["length_consistent"]  # Should be flagged
