from __future__ import annotations

"""Comprehensive edge case tests for BioSeqFlow."""


from typing import TYPE_CHECKING

import pytest
from pydantic import ValidationError

from bioseqflow.core.config import Config
from bioseqflow.preprocessing.deduplication import DuplicateRemover, UMIDeduplicator
from bioseqflow.preprocessing.filtering import LengthFilter, QualityFilter
from bioseqflow.preprocessing.trimming import AdapterTrimmer, QualityTrimmer
from bioseqflow.utils.alignment import find_adapter_fuzzy, smith_waterman
from bioseqflow.utils.io import FastqRecord, read_fastq, write_fastq

if TYPE_CHECKING:
    from pathlib import Path


class TestFastqEdgeCases:
    """Test FASTQ I/O edge cases."""

    def test_empty_fastq_file(self, tmp_path: Path):
        """Test reading empty FASTQ file."""
        empty_file = tmp_path / "empty.fastq"
        empty_file.write_text("")

        records = list(read_fastq(empty_file))
        assert len(records) == 0

    def test_single_read(self, tmp_path: Path):
        """Test file with single read."""
        single_file = tmp_path / "single.fastq"
        records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
        ]
        write_fastq(records, single_file)

        read_records = list(read_fastq(single_file))
        assert len(read_records) == 1

    def test_very_long_read(self, tmp_path: Path):
        """Test with very long read (>10kb)."""
        long_file = tmp_path / "long.fastq"
        sequence = "A" * 15000
        quality = "I" * 15000

        records = [
            FastqRecord("@long_read\n", sequence + "\n", "+\n", quality + "\n"),
        ]
        write_fastq(records, long_file)

        read_records = list(read_fastq(long_file))
        assert len(read_records) == 1
        assert read_records[0].length == 15000

    def test_read_with_special_characters(self, tmp_path: Path):
        """Test read with special characters in header."""
        special_file = tmp_path / "special.fastq"
        records = [
            FastqRecord(
                "@read:1:2:3 special/chars\n",
                "ATCG\n",
                "+\n",
                "IIII\n",
            ),
        ]
        write_fastq(records, special_file)

        read_records = list(read_fastq(special_file))
        assert len(read_records) == 1

    def test_all_n_sequence(self, tmp_path: Path):
        """Test read with all N bases."""
        n_file = tmp_path / "all_n.fastq"
        records = [
            FastqRecord("@read_n\n", "NNNN\n", "+\n", "####\n"),
        ]
        write_fastq(records, n_file)

        read_records = list(read_fastq(n_file))
        assert len(read_records) == 1
        assert read_records[0].sequence.strip() == "NNNN"

    def test_minimum_length_read(self, tmp_path: Path):
        """Test read with minimum length (1bp)."""
        min_file = tmp_path / "min.fastq"
        records = [
            FastqRecord("@min\n", "A\n", "+\n", "I\n"),
        ]
        write_fastq(records, min_file)

        read_records = list(read_fastq(min_file))
        assert len(read_records) == 1
        assert read_records[0].length == 1


class TestQualityFilterEdgeCases:
    """Test quality filtering edge cases."""

    def test_all_low_quality_file(self, tmp_path: Path):
        """Test file where all reads fail quality."""
        low_q_file = tmp_path / "low_quality.fastq"
        records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "####\n"),  # All Q2
            FastqRecord("@read2\n", "GCTA\n", "+\n", "!!!!  \n"),  # All Q0-1
        ]
        write_fastq(records, low_q_file)

        output = tmp_path / "filtered.fastq"
        filter_obj = QualityFilter()
        stats = filter_obj.filter(low_q_file, output, min_quality=20)

        assert stats["passed_reads"] == 0
        assert stats["failed_quality"] == 2

    def test_zero_length_after_filtering(self, tmp_path: Path):
        """Test read that becomes zero length after filtering."""
        input_file = tmp_path / "input.fastq"
        records = [
            FastqRecord("@short\n", "AT\n", "+\n", "II\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "filtered.fastq"
        filter_obj = LengthFilter()
        stats = filter_obj.filter(input_file, output, min_length=10)

        assert stats["too_short"] == 1
        assert stats["passed_reads"] == 0

    def test_max_n_content_boundary(self, tmp_path: Path):
        """Test N content at exact boundary."""
        input_file = tmp_path / "n_content.fastq"
        # 20% N content: 10 N's in 50bp read
        seq = "ATCGATCGATNNNNNNNNNNGATCGATCGATCGATCGATCGATCGATCGA"  # 50bp
        assert len(seq) == 50, f"Sequence length is {len(seq)}, not 50"
        qual = "I" * len(seq)
        records = [
            FastqRecord("@n_test\n", seq + "\n", "+\n", qual + "\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "filtered.fastq"
        filter_obj = QualityFilter()

        # Should pass with 20% threshold (exactly at boundary)
        stats = filter_obj.filter(input_file, output, min_quality=0, min_length=10, max_n_content=20.0)
        assert stats["passed_reads"] == 1

        # Should pass with 21% threshold
        output2 = tmp_path / "filtered2.fastq"
        stats2 = filter_obj.filter(input_file, output2, min_quality=0, min_length=10, max_n_content=21.0)
        assert stats2["passed_reads"] == 1

        # Should fail with 19% threshold
        output3 = tmp_path / "filtered3.fastq"
        stats3 = filter_obj.filter(input_file, output3, min_quality=0, min_length=10, max_n_content=19.0)
        assert stats3["passed_reads"] == 0

    def test_quality_score_boundaries(self, tmp_path: Path):
        """Test quality scores at exact thresholds."""
        input_file = tmp_path / "quality_boundary.fastq"
        # Use a longer read to ensure it passes length filter
        # Mean quality exactly 20 (ASCII '5' = 53, 53-33 = 20)
        seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"  # 60bp
        qual = "5" * 60  # All Q20
        records = [
            FastqRecord("@q20\n", seq + "\n", "+\n", qual + "\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "filtered.fastq"
        filter_obj = QualityFilter()

        # Should pass with threshold 20 (>= comparison)
        stats = filter_obj.filter(input_file, output, min_quality=20, min_length=50)
        assert stats["passed_reads"] == 1

        # Should fail with threshold 21 (mean is 20)
        output2 = tmp_path / "filtered2.fastq"
        stats2 = filter_obj.filter(input_file, output2, min_quality=21, min_length=50)
        assert stats2["passed_reads"] == 0


class TestDeduplicationEdgeCases:
    """Test deduplication edge cases."""

    def test_all_unique_reads(self, tmp_path: Path):
        """Test file with no duplicates."""
        input_file = tmp_path / "unique.fastq"
        records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
            FastqRecord("@read2\n", "GCTA\n", "+\n", "IIII\n"),
            FastqRecord("@read3\n", "TTAA\n", "+\n", "IIII\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "dedup.fastq"
        dedup = DuplicateRemover()
        stats = dedup.remove_duplicates(input_file, output)

        assert stats["unique_reads"] == 3
        assert stats["duplicates"] == 0

    def test_all_duplicate_reads(self, tmp_path: Path):
        """Test file where all reads are identical."""
        input_file = tmp_path / "all_dup.fastq"
        records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),
            FastqRecord("@read2\n", "ATCG\n", "+\n", "IIII\n"),
            FastqRecord("@read3\n", "ATCG\n", "+\n", "IIII\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "dedup.fastq"
        dedup = DuplicateRemover()
        stats = dedup.remove_duplicates(input_file, output)

        assert stats["unique_reads"] == 1
        assert stats["duplicates"] == 2

    def test_hash_collision_handling(self, tmp_path: Path):
        """Test that hash method handles many duplicates."""
        input_file = tmp_path / "many_dup.fastq"
        # Create 1000 copies of same sequence
        records = [
            FastqRecord(f"@read{i}\n", "ATCGATCG\n", "+\n", "IIIIIIII\n")
            for i in range(1000)
        ]
        write_fastq(records, input_file)

        output = tmp_path / "dedup.fastq"
        dedup = DuplicateRemover()
        stats = dedup.remove_duplicates(input_file, output, use_hash=True)

        assert stats["unique_reads"] == 1
        assert stats["duplicates"] == 999

    def test_prefix_dedup_short_reads(self, tmp_path: Path):
        """Test prefix deduplication with reads shorter than prefix."""
        input_file = tmp_path / "short.fastq"
        records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),  # 4bp < 30bp prefix
            FastqRecord("@read2\n", "ATCG\n", "+\n", "IIII\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "dedup.fastq"
        dedup = DuplicateRemover()
        stats = dedup.remove_duplicates(input_file, output, method="prefix")

        assert stats["unique_reads"] == 1


class TestUMIEdgeCases:
    """Test UMI deduplication edge cases."""

    def test_umi_longer_than_read(self, tmp_path: Path):
        """Test UMI length longer than read."""
        input_file = tmp_path / "short_read.fastq"
        records = [
            FastqRecord("@read1\n", "ATCG\n", "+\n", "IIII\n"),  # 4bp read
        ]
        write_fastq(records, input_file)

        output = tmp_path / "dedup.fastq"
        dedup = UMIDeduplicator()

        # Should handle gracefully (UMI length 8 > read length 4)
        stats = dedup.deduplicate(input_file, output, umi_length=8)
        # All sequence would be UMI, nothing left
        assert stats["total_reads"] == 1

    def test_zero_edit_distance(self, tmp_path: Path):
        """Test UMI error correction with zero edit distance (exact match)."""
        input_file = tmp_path / "umi.fastq"
        records = [
            FastqRecord("@read1\n", "ATCGATCGATCG\n", "+\n", "IIIIIIIIIIII\n"),
            FastqRecord("@read2\n", "ATCGATCGATCG\n", "+\n", "IIIIIIIIIIII\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "dedup.fastq"
        dedup = UMIDeduplicator()
        stats = dedup.deduplicate(
            input_file, output, umi_length=4, error_correction=True, max_edit_distance=0
        )

        assert stats["unique_umis"] == 1

    def test_all_different_umis(self, tmp_path: Path):
        """Test when all UMIs are unique."""
        input_file = tmp_path / "unique_umis.fastq"
        records = [
            FastqRecord("@read1\n", "AAAATCGATCG\n", "+\n", "IIIIIIIIIII\n"),
            FastqRecord("@read2\n", "TTTATCGATCG\n", "+\n", "IIIIIIIIIII\n"),
            FastqRecord("@read3\n", "GGGGATCGATCG\n", "+\n", "IIIIIIIIIIII\n"),  # 12 I's
        ]
        write_fastq(records, input_file)

        output = tmp_path / "dedup.fastq"
        dedup = UMIDeduplicator()
        stats = dedup.deduplicate(input_file, output, umi_length=4)

        assert stats["unique_umis"] == 3
        assert stats["duplicates"] == 0

    def test_high_error_rate_correction(self, tmp_path: Path):
        """Test UMI correction with very similar UMIs."""
        input_file = tmp_path / "similar_umis.fastq"
        # UMIs differ by 1 base
        records = [
            FastqRecord("@read1\n", "AAAATCGATCG\n", "+\n", "IIIIIIIIIII\n"),
            FastqRecord("@read2\n", "AAAATCGATCG\n", "+\n", "IIIIIIIIIII\n"),
            FastqRecord("@read3\n", "AAAATCGATCG\n", "+\n", "IIIIIIIIIII\n"),
            FastqRecord("@read4\n", "AAATTCGATCG\n", "+\n", "IIIIIIIIIII\n"),  # 1 diff
        ]
        write_fastq(records, input_file)

        output = tmp_path / "dedup.fastq"
        dedup = UMIDeduplicator()
        stats = dedup.deduplicate(
            input_file, output, umi_length=4, error_correction=True, max_edit_distance=1
        )

        # Should merge the similar UMI
        assert stats["umis_merged"] > 0


class TestAlignmentEdgeCases:
    """Test alignment edge cases."""

    def test_empty_sequences(self):
        """Test alignment with empty sequences."""
        result = smith_waterman("", "ATCG")
        assert result.score == 0

        result = smith_waterman("ATCG", "")
        assert result.score == 0

        result = smith_waterman("", "")
        assert result.score == 0

    def test_identical_sequences(self):
        """Test alignment of identical sequences."""
        seq = "ATCGATCGATCG"
        result = smith_waterman(seq, seq)

        assert result.score > 0
        assert result.aligned_query == seq
        assert result.aligned_target == seq

    def test_no_similarity(self):
        """Test sequences with no similarity."""
        result = smith_waterman("AAAA", "TTTT")

        # Should have very low or zero score
        assert result.score >= 0

    def test_adapter_longer_than_read(self):
        """Test finding adapter longer than read."""
        adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # Long adapter (35bp)
        read = "TTTTTTTT"  # Short read with no match

        result = find_adapter_fuzzy(read, adapter, min_overlap=5)
        # Should not find it (no significant match)
        assert result is None

    def test_very_long_sequences(self):
        """Test alignment with very long sequences (performance test)."""
        long_query = "A" * 1000
        long_target = "A" * 1000

        # Should complete without timeout/memory issues
        result = smith_waterman(long_query, long_target)
        assert result.score > 0

    def test_single_base_sequences(self):
        """Test alignment with single base sequences."""
        result = smith_waterman("A", "A")
        assert result.score > 0

        result = smith_waterman("A", "T")
        assert result.score >= 0


class TestConfigEdgeCases:
    """Test configuration edge cases."""

    def test_invalid_organism(self):
        """Test config with invalid organism."""
        with pytest.raises(ValueError, match="Unknown organism"):
            Config(organism="invalid_species")

    def test_quality_score_boundaries(self):
        """Test quality score at boundaries."""
        # Valid boundaries
        config = Config(min_quality=0)
        assert config.min_quality == 0

        config = Config(min_quality=93)
        assert config.min_quality == 93

        # Invalid boundaries should raise validation error
        with pytest.raises(ValidationError):
            Config(min_quality=-1)

        with pytest.raises(ValidationError):
            Config(min_quality=94)

    def test_gc_content_boundaries(self):
        """Test GC content at exact boundaries."""
        config = Config(gc_content_min=0.0, gc_content_max=1.0)
        assert config.gc_content_min == 0.0
        assert config.gc_content_max == 1.0

        # Invalid ranges
        with pytest.raises(ValidationError):
            Config(gc_content_min=-0.1)

        with pytest.raises(ValidationError):
            Config(gc_content_max=1.1)

    def test_organism_overrides_gc(self):
        """Test that organism profile overrides GC settings."""
        config = Config(
            organism="malaria",
            gc_content_min=0.4,  # Will be overridden
            gc_content_max=0.6,  # Will be overridden
        )

        # Should use malaria profile (0.15, 0.25)
        assert config.gc_content_min == 0.15
        assert config.gc_content_max == 0.25


class TestAdapterTrimmingEdgeCases:
    """Test adapter trimming edge cases."""

    def test_adapter_at_read_start(self, tmp_path: Path):
        """Test adapter at the very start of read."""
        input_file = tmp_path / "adapter_start.fastq"
        adapter = "AGATCGGAAGAGC"
        records = [
            FastqRecord("@read1\n", adapter + "ATCGATCG\n", "+\n", "I" * 21 + "\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "trimmed.fastq"
        trimmer = AdapterTrimmer()
        stats = trimmer.trim(input_file, output, adapter, use_cutadapt=False)

        # Should find adapter at position 0
        assert stats["trimmed_reads"] >= 1

    def test_multiple_adapter_occurrences(self, tmp_path: Path):
        """Test read with adapter appearing multiple times."""
        input_file = tmp_path / "multi_adapter.fastq"
        adapter = "AGATCGGAAGAGC"
        records = [
            FastqRecord(
                "@read1\n",
                "ATCG" + adapter + "GCTA" + adapter + "\n",
                "+\n",
                "I" * 34 + "\n",
            ),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "trimmed.fastq"
        trimmer = AdapterTrimmer()
        stats = trimmer.trim(input_file, output, adapter, use_cutadapt=False)

        # Should find first occurrence
        assert stats["trimmed_reads"] == 1

    def test_adapter_equals_read(self, tmp_path: Path):
        """Test when adapter is entire read."""
        input_file = tmp_path / "adapter_only.fastq"
        adapter = "AGATCGGAAGAGC"
        records = [
            FastqRecord("@read1\n", adapter + "\n", "+\n", "I" * 13 + "\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "trimmed.fastq"
        trimmer = AdapterTrimmer()
        stats = trimmer.trim(input_file, output, adapter, use_cutadapt=False)

        # Read should be completely trimmed
        assert stats["trimmed_reads"] == 1


class TestQualityTrimmingEdgeCases:
    """Test quality trimming edge cases."""

    def test_all_high_quality(self, tmp_path: Path):
        """Test read with all high quality bases."""
        input_file = tmp_path / "high_q.fastq"
        records = [
            FastqRecord("@read1\n", "ATCGATCG\n", "+\n", "IIIIIIII\n"),  # All Q40
        ]
        write_fastq(records, input_file)

        output = tmp_path / "trimmed.fastq"
        trimmer = QualityTrimmer()
        stats = trimmer.trim(input_file, output, min_quality=20)

        # Nothing should be trimmed
        assert stats["trimmed_reads"] == 0

    def test_alternating_quality(self, tmp_path: Path):
        """Test read with alternating high/low quality."""
        input_file = tmp_path / "alt_q.fastq"
        records = [
            # Alternating Q40 (I) and Q2 (#)
            FastqRecord("@read1\n", "ATCGATCG\n", "+\n", "I#I#I#I#\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "trimmed.fastq"
        trimmer = QualityTrimmer()
        stats = trimmer.trim(input_file, output, min_quality=20, trim_5prime=True)

        # Should trim something
        assert stats["total_reads"] == 1

    def test_single_base_read(self, tmp_path: Path):
        """Test trimming single base read."""
        input_file = tmp_path / "single_base.fastq"
        records = [
            FastqRecord("@read1\n", "A\n", "+\n", "I\n"),
        ]
        write_fastq(records, input_file)

        output = tmp_path / "trimmed.fastq"
        trimmer = QualityTrimmer()
        stats = trimmer.trim(input_file, output, min_quality=20)

        # Should keep the single base if quality is good
        assert stats["total_reads"] == 1
