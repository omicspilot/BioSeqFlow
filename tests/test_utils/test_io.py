"""Tests for I/O utilities."""

from pathlib import Path

import pytest

from bioseqflow.utils.io import (
    FastqRecord,
    count_reads,
    get_read_length_distribution,
    open_file,
    parse_sample_sheet,
    read_fastq,
    write_fastq,
)


class TestFastqRecord:
    """Tests for FastqRecord class."""

    def test_fastq_record_creation(self) -> None:
        """Test creating a FASTQ record."""
        record = FastqRecord(
            "@SEQ_ID", "GATTTGGGGTTCAAAGCAGTATCG", "+", "!''*((((***+))%%%++)(%%"
        )

        assert record.header == "@SEQ_ID"
        assert record.sequence == "GATTTGGGGTTCAAAGCAGTATCG"
        assert record.quality == "!''*((((***+))%%%++)(%%"

    def test_fastq_record_name(self) -> None:
        """Test extracting sequence name."""
        record = FastqRecord("@SEQ_ID_1", "GATC", "+", "!!!!")

        assert record.name == "SEQ_ID_1"

    def test_fastq_record_length(self) -> None:
        """Test getting sequence length."""
        record = FastqRecord("@SEQ_ID", "GATTTGGGGTTCAAAGCAGTATCG", "+", "!" * 24)

        assert record.length == 24

    def test_fastq_record_mean_quality(self) -> None:
        """Test calculating mean quality."""
        # Quality string with known values
        # "!" = Phred 0, "I" = Phred 40
        record = FastqRecord("@SEQ_ID", "GATC", "+", "!!!!")

        mean_qual = record.mean_quality()
        assert mean_qual == 0.0

    def test_fastq_record_str(self) -> None:
        """Test string representation."""
        record = FastqRecord("@SEQ_ID", "GATC", "+", "!!!!")

        string_repr = str(record)
        assert string_repr == "@SEQ_ID\nGATC\n+\n!!!!\n"


class TestFileIO:
    """Tests for file I/O operations."""

    def test_open_file_regular(self, temp_dir: Path) -> None:
        """Test opening regular file."""
        test_file = temp_dir / "test.txt"
        test_file.write_text("test content")

        with open_file(test_file, "r") as f:
            content = f.read()

        assert content == "test content"

    def test_open_file_gzip(self, temp_dir: Path) -> None:
        """Test opening gzip file."""
        test_file = temp_dir / "test.txt.gz"

        with open_file(test_file, "w") as f:
            f.write("compressed content")

        with open_file(test_file, "r") as f:
            content = f.read()

        assert content == "compressed content"

    def test_open_file_not_found(self, temp_dir: Path) -> None:
        """Test opening non-existent file."""
        with pytest.raises(FileNotFoundError), open_file(temp_dir / "nonexistent.txt", "r") as _f:
            pass

    def test_read_fastq(self, sample_fastq_file: Path) -> None:
        """Test reading FASTQ file."""
        records = list(read_fastq(sample_fastq_file))

        assert len(records) == 2
        assert records[0].header.startswith("@SEQ_ID")
        assert len(records[0].sequence) > 0
        assert len(records[0].quality) == len(records[0].sequence)

    def test_write_fastq(self, temp_dir: Path) -> None:
        """Test writing FASTQ file."""
        records = [
            FastqRecord("@SEQ1", "GATC", "+", "!!!!"),
            FastqRecord("@SEQ2", "ATCG", "+", "IIII"),
        ]

        output_file = temp_dir / "output.fastq"
        count = write_fastq(records, output_file)

        assert count == 2
        assert output_file.exists()

        # Read back and verify
        read_records = list(read_fastq(output_file))
        assert len(read_records) == 2
        assert read_records[0].sequence == "GATC"

    def test_write_fastq_compressed(self, temp_dir: Path) -> None:
        """Test writing compressed FASTQ."""
        records = [FastqRecord("@SEQ1", "GATC", "+", "!!!!")]

        output_file = temp_dir / "output.fastq"
        write_fastq(records, output_file, compress=True)

        # Should create .gz file
        assert (temp_dir / "output.fastq.gz").exists()

    def test_count_reads(self, sample_fastq_file: Path) -> None:
        """Test counting reads."""
        count = count_reads(sample_fastq_file)

        assert count == 2

    def test_get_read_length_distribution(self, sample_fastq_file: Path) -> None:
        """Test getting read length distribution."""
        dist = get_read_length_distribution(sample_fastq_file)

        assert isinstance(dist, dict)
        assert len(dist) > 0
        assert all(isinstance(k, int) for k in dist)
        assert all(isinstance(v, int) for v in dist.values())


class TestSampleSheet:
    """Tests for sample sheet parsing."""

    def test_parse_sample_sheet(self, temp_dir: Path) -> None:
        """Test parsing sample sheet."""
        csv_file = temp_dir / "samples.csv"
        csv_content = """sample_id,fastq_r1,fastq_r2,adapter
sample1,s1_R1.fastq.gz,s1_R2.fastq.gz,AGATCGGAAGAG
sample2,s2_R1.fastq.gz,s2_R2.fastq.gz,AGATCGGAAGAG
"""
        csv_file.write_text(csv_content)

        samples = parse_sample_sheet(csv_file)

        assert len(samples) == 2
        assert samples[0]["sample_id"] == "sample1"
        assert samples[0]["fastq_r1"] == "s1_R1.fastq.gz"
        assert samples[1]["adapter"] == "AGATCGGAAGAG"

    def test_parse_sample_sheet_missing_file(self, temp_dir: Path) -> None:
        """Test parsing non-existent sample sheet."""
        with pytest.raises(FileNotFoundError):
            parse_sample_sheet(temp_dir / "nonexistent.csv")

    def test_parse_sample_sheet_missing_columns(self, temp_dir: Path) -> None:
        """Test parsing sample sheet with missing required columns."""
        csv_file = temp_dir / "bad_samples.csv"
        csv_content = """sample_id
sample1
"""
        csv_file.write_text(csv_content)

        with pytest.raises(ValueError, match="Missing required columns"):
            parse_sample_sheet(csv_file)
