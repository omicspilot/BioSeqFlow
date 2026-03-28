from __future__ import annotations

"""Tests for QC parsers."""


from typing import TYPE_CHECKING

import pytest

from bioseqflow.qc.parsers import FastQCParser, MultiQCParser

if TYPE_CHECKING:
    from pathlib import Path


class TestFastQCParser:
    """Test FastQC parser."""

    def test_parse_basic_fastqc_file(self, tmp_path: Path):
        """Test parsing basic FastQC data file."""
        fastqc_file = tmp_path / "fastqc_data.txt"

        # Create sample FastQC output
        content = """##FastQC	0.11.9
>>Basic Statistics	pass
#Measure	Value
Filename	sample.fastq
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	1000
Sequences flagged as poor quality	0
Sequence length	50
%GC	45
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	35.5	36	35	37	34	38
2	35.3	36	35	37	33	38
>>END_MODULE
>>Per sequence quality scores	warn
#Quality	Count
20	10
30	50
40	940
>>END_MODULE
"""
        fastqc_file.write_text(content)

        parser = FastQCParser()
        results = parser.parse(fastqc_file)

        # Check structure
        assert "basic_statistics" in results
        assert "modules" in results
        assert "summary" in results

        # Check summary statuses
        assert results["summary"]["Basic Statistics"] == "pass"
        assert results["summary"]["Per base sequence quality"] == "pass"
        assert results["summary"]["Per sequence quality scores"] == "warn"

        # Check Basic Statistics module
        basic_stats = results["modules"]["Basic Statistics"]
        assert basic_stats["Filename"] == "sample.fastq"
        assert basic_stats["Total Sequences"] == 1000
        assert basic_stats["%GC"] == 45
        assert basic_stats["Encoding"] == "Sanger / Illumina 1.9"

        # Check other modules have data
        assert "Per base sequence quality" in results["modules"]
        assert "Per sequence quality scores" in results["modules"]

    def test_parse_file_not_found(self, tmp_path: Path):
        """Test error when FastQC file doesn't exist."""
        parser = FastQCParser()

        with pytest.raises(FileNotFoundError, match="FastQC data file not found"):
            parser.parse(tmp_path / "nonexistent.txt")

    def test_parse_with_comments(self, tmp_path: Path):
        """Test parsing file with comment lines."""
        fastqc_file = tmp_path / "fastqc_data.txt"

        content = """##FastQC	0.11.9
# This is a comment
>>Basic Statistics	pass
#Measure	Value
Filename	test.fastq
Total Sequences	500
>>END_MODULE
"""
        fastqc_file.write_text(content)

        parser = FastQCParser()
        results = parser.parse(fastqc_file)

        assert results["modules"]["Basic Statistics"]["Total Sequences"] == 500

    def test_parse_basic_statistics_numeric_conversion(self, tmp_path: Path):
        """Test that numeric values are converted to int."""
        fastqc_file = tmp_path / "fastqc_data.txt"

        content = """##FastQC	0.11.9
>>Basic Statistics	pass
#Measure	Value
Total Sequences	1234
%GC	42
Filename	test.fastq
>>END_MODULE
"""
        fastqc_file.write_text(content)

        parser = FastQCParser()
        results = parser.parse(fastqc_file)

        stats = results["modules"]["Basic Statistics"]
        assert isinstance(stats["Total Sequences"], int)
        assert isinstance(stats["%GC"], int)
        assert isinstance(stats["Filename"], str)

    def test_parse_module_without_summary(self, tmp_path: Path):
        """Test parsing module that has no summary status."""
        fastqc_file = tmp_path / "fastqc_data.txt"

        content = """##FastQC	0.11.9
>>Custom Module
#Header
data line 1
data line 2
>>END_MODULE
"""
        fastqc_file.write_text(content)

        parser = FastQCParser()
        results = parser.parse(fastqc_file)

        assert "Custom Module" in results["modules"]
        assert results["modules"]["Custom Module"]["raw_data"] == ["data line 1", "data line 2"]

    def test_parse_empty_module(self, tmp_path: Path):
        """Test parsing module with no data lines."""
        fastqc_file = tmp_path / "fastqc_data.txt"

        content = """##FastQC	0.11.9
>>Empty Module	pass
>>END_MODULE
>>Basic Statistics	pass
#Measure	Value
Total Sequences	100
>>END_MODULE
"""
        fastqc_file.write_text(content)

        parser = FastQCParser()
        results = parser.parse(fastqc_file)

        # Empty module should not be in results
        assert "Empty Module" not in results["modules"]
        # But summary should still be there
        assert results["summary"]["Empty Module"] == "pass"


class TestMultiQCParser:
    """Test MultiQC parser."""

    def test_parse_multiqc_directory(self, tmp_path: Path):
        """Test parsing MultiQC data directory."""
        multiqc_dir = tmp_path / "multiqc_data"
        multiqc_dir.mkdir()

        # Create general stats file
        general_stats = multiqc_dir / "multiqc_general_stats.txt"
        general_stats.write_text(
            """Sample	Total Sequences	%GC	Avg Quality
sample1	1000	45	35.5
sample2	2000	42	36.2
"""
        )

        # Create FastQC data file
        fastqc_data = multiqc_dir / "multiqc_fastqc.txt"
        fastqc_data.write_text(
            """Sample	total_sequences	percent_gc	avg_sequence_length
sample1	1000	45.0	150
sample2	2000	42.0	150
"""
        )

        parser = MultiQCParser()
        results = parser.parse(multiqc_dir)

        # Check general stats
        assert "sample1" in results["general_stats"]
        assert results["general_stats"]["sample1"]["Total Sequences"] == "1000"
        assert results["general_stats"]["sample1"]["%GC"] == "45"

        assert "sample2" in results["general_stats"]
        assert results["general_stats"]["sample2"]["Total Sequences"] == "2000"

        # Check FastQC data
        assert "sample1" in results["fastqc_data"]
        assert results["fastqc_data"]["sample1"]["total_sequences"] == "1000"
        assert results["fastqc_data"]["sample2"]["percent_gc"] == "42.0"

    def test_parse_directory_not_found(self, tmp_path: Path):
        """Test error when MultiQC directory doesn't exist."""
        parser = MultiQCParser()

        with pytest.raises(FileNotFoundError, match="MultiQC data directory not found"):
            parser.parse(tmp_path / "nonexistent_dir")

    def test_parse_directory_without_files(self, tmp_path: Path):
        """Test parsing directory with no data files."""
        multiqc_dir = tmp_path / "multiqc_data"
        multiqc_dir.mkdir()

        parser = MultiQCParser()
        results = parser.parse(multiqc_dir)

        # Should return empty dicts
        assert results["general_stats"] == {}
        assert results["fastqc_data"] == {}

    def test_parse_only_general_stats(self, tmp_path: Path):
        """Test parsing with only general stats file."""
        multiqc_dir = tmp_path / "multiqc_data"
        multiqc_dir.mkdir()

        general_stats = multiqc_dir / "multiqc_general_stats.txt"
        general_stats.write_text(
            """Sample	Total Sequences
sample1	1000
"""
        )

        parser = MultiQCParser()
        results = parser.parse(multiqc_dir)

        assert len(results["general_stats"]) == 1
        assert results["fastqc_data"] == {}

    def test_parse_only_fastqc_data(self, tmp_path: Path):
        """Test parsing with only FastQC data file."""
        multiqc_dir = tmp_path / "multiqc_data"
        multiqc_dir.mkdir()

        fastqc_data = multiqc_dir / "multiqc_fastqc.txt"
        fastqc_data.write_text(
            """Sample	total_sequences
sample1	1000
sample2	2000
"""
        )

        parser = MultiQCParser()
        results = parser.parse(multiqc_dir)

        assert results["general_stats"] == {}
        assert len(results["fastqc_data"]) == 2

    def test_parse_multiple_samples(self, tmp_path: Path):
        """Test parsing with multiple samples."""
        multiqc_dir = tmp_path / "multiqc_data"
        multiqc_dir.mkdir()

        general_stats = multiqc_dir / "multiqc_general_stats.txt"
        general_stats.write_text(
            """Sample	Reads	GC
sample1	100	45
sample2	200	42
sample3	300	48
sample4	400	50
"""
        )

        parser = MultiQCParser()
        results = parser.parse(multiqc_dir)

        assert len(results["general_stats"]) == 4
        assert results["general_stats"]["sample1"]["Reads"] == "100"
        assert results["general_stats"]["sample4"]["GC"] == "50"
