from __future__ import annotations

"""Tests for fuzzy adapter trimming."""

from pathlib import Path

import pytest

from bioseqflow.preprocessing.trimming import AdapterTrimmer
from bioseqflow.utils.io import FastqRecord, write_fastq


@pytest.fixture
def temp_fastq_fuzzy(tmp_path: Path) -> Path:
    """Create temporary FASTQ file with adapter sequences."""
    fastq_file = tmp_path / "test_fuzzy.fastq"

    # Create test records with adapters (some with errors)
    records = [
        # Perfect adapter match
        FastqRecord(
            "@read1\n",
            "ATCGATCGAGATCGGAAGAGCGTCGT\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
        # Adapter with 1 mismatch (G->T at position 7)
        FastqRecord(
            "@read2\n",
            "ATCGATCGAGATCGTAAGAGCGTCGT\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
        # Partial adapter at end
        FastqRecord(
            "@read3\n",
            "ATCGATCGATCGATCGAGATCGGAAG\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
        # No adapter
        FastqRecord(
            "@read4\n",
            "ATCGATCGATCGATCGATCGATCG\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
    ]

    write_fastq(records, fastq_file)
    return fastq_file


class TestFuzzyAdapterTrimming:
    """Test fuzzy adapter trimming functionality."""

    def test_exact_match_fuzzy_mode(self, temp_fastq_fuzzy: Path, tmp_path: Path):
        """Test that fuzzy matching still finds exact matches."""
        output_file = tmp_path / "trimmed_exact.fastq"
        adapter = "AGATCGGAAGAGC"

        trimmer = AdapterTrimmer()
        stats = trimmer.trim(
            temp_fastq_fuzzy,
            output_file,
            adapter,
            use_cutadapt=False,
            fuzzy_match=True,
        )

        # Should trim at least the perfect match
        assert stats["trimmed_reads"] >= 1
        assert stats["total_bp_removed"] > 0

    def test_mismatch_detection(self, temp_fastq_fuzzy: Path, tmp_path: Path):
        """Test that fuzzy matching detects adapters with mismatches."""
        output_file = tmp_path / "trimmed_mismatch.fastq"
        adapter = "AGATCGGAAGAGC"

        # Fuzzy matching should find more adapters than exact
        trimmer = AdapterTrimmer()

        # Exact matching
        exact_output = tmp_path / "exact.fastq"
        exact_stats = trimmer.trim(
            temp_fastq_fuzzy,
            exact_output,
            adapter,
            use_cutadapt=False,
            fuzzy_match=False,
        )

        # Fuzzy matching
        fuzzy_stats = trimmer.trim(
            temp_fastq_fuzzy,
            output_file,
            adapter,
            use_cutadapt=False,
            fuzzy_match=True,
            max_error_rate=0.15,
        )

        # Fuzzy should find at least as many as exact (likely more)
        assert fuzzy_stats["trimmed_reads"] >= exact_stats["trimmed_reads"]

    def test_error_rate_threshold(self, temp_fastq_fuzzy: Path, tmp_path: Path):
        """Test different error rate thresholds."""
        adapter = "AGATCGGAAGAGC"
        trimmer = AdapterTrimmer()

        # Strict threshold
        strict_output = tmp_path / "strict.fastq"
        strict_stats = trimmer.trim(
            temp_fastq_fuzzy,
            strict_output,
            adapter,
            use_cutadapt=False,
            fuzzy_match=True,
            max_error_rate=0.05,
        )

        # Lenient threshold
        lenient_output = tmp_path / "lenient.fastq"
        lenient_stats = trimmer.trim(
            temp_fastq_fuzzy,
            lenient_output,
            adapter,
            use_cutadapt=False,
            fuzzy_match=True,
            max_error_rate=0.20,
        )

        # Lenient should trim at least as many as strict
        assert lenient_stats["trimmed_reads"] >= strict_stats["trimmed_reads"]

    def test_partial_adapter_fuzzy(self, tmp_path: Path):
        """Test fuzzy matching with partial adapter."""
        # Create file with partial adapter at end
        fastq_file = tmp_path / "partial.fastq"
        records = [
            FastqRecord(
                "@read_partial\n",
                "ATCGATCGATCGATCGAGATCGGAA\n",  # Only first 8bp of adapter
                "+\n",
                "IIIIIIIIIIIIIIIIIIIIIIIII\n",
            ),
        ]
        write_fastq(records, fastq_file)

        output_file = tmp_path / "trimmed_partial.fastq"
        adapter = "AGATCGGAAGAGC"

        trimmer = AdapterTrimmer()
        trimmer.min_overlap = 6  # Lower threshold for partial match
        stats = trimmer.trim(
            fastq_file,
            output_file,
            adapter,
            use_cutadapt=False,
            fuzzy_match=True,
        )

        # Should find and trim partial adapter
        assert stats["trimmed_reads"] >= 1

    def test_no_false_positives(self, tmp_path: Path):
        """Test that fuzzy matching doesn't create false positives."""
        # Create file with no real adapter
        fastq_file = tmp_path / "no_adapter.fastq"
        records = [
            FastqRecord(
                "@read_clean\n",
                "ATCGATCGATCGATCGATCGATCG\n",
                "+\n",
                "IIIIIIIIIIIIIIIIIIIIIIII\n",
            ),
        ]
        write_fastq(records, fastq_file)

        output_file = tmp_path / "trimmed_clean.fastq"
        adapter = "AGATCGGAAGAGC"

        trimmer = AdapterTrimmer()
        stats = trimmer.trim(
            fastq_file,
            output_file,
            adapter,
            use_cutadapt=False,
            fuzzy_match=True,
            max_error_rate=0.15,
        )

        # Should not trim reads without adapter
        assert stats["trimmed_reads"] == 0
