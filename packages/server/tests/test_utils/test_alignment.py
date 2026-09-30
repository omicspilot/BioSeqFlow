from __future__ import annotations

"""Tests for sequence alignment utilities."""


from bioseqflow_server.toolkit.utils.alignment import (
    calculate_match_score,
    find_adapter_fuzzy,
    smith_waterman,
)


class TestSmithWaterman:
    """Test Smith-Waterman alignment algorithm."""

    def test_perfect_match(self):
        """Test perfect alignment match."""
        query = "ATCG"
        target = "GGATCGCC"

        result = smith_waterman(query, target)

        assert result.score > 0
        assert result.aligned_query == "ATCG"
        assert result.aligned_target == "ATCG"
        assert result.start == 2
        assert result.end == 6

    def test_mismatch_alignment(self):
        """Test alignment with mismatches."""
        query = "ATCG"
        target = "GGATTGCC"  # TC -> TT (one mismatch)

        result = smith_waterman(query, target)

        assert result.score > 0
        # Should still find the best local alignment
        assert len(result.aligned_query) > 0

    def test_gap_alignment(self):
        """Test alignment with gaps."""
        query = "ATCG"
        target = "GGATGCC"  # Missing C

        result = smith_waterman(query, target)

        assert result.score > 0
        assert len(result.aligned_query) > 0

    def test_no_match(self):
        """Test with no significant match."""
        query = "ATCG"
        target = "GGGG"

        result = smith_waterman(query, target)

        # Should have low or zero score
        assert result.score >= 0

    def test_adapter_in_middle(self):
        """Test finding adapter in middle of sequence."""
        adapter = "AGATCGGAAGAGC"
        read = "ATCGATCGAGATCGGAAGAGCGTCGT"

        result = smith_waterman(adapter, read)

        assert result.score > 0
        assert result.start == 8
        # Should align the adapter sequence


class TestFindAdapterFuzzy:
    """Test fuzzy adapter finding."""

    def test_exact_adapter_match(self):
        """Test finding exact adapter match."""
        sequence = "ATCGATCGAGATCGGAAGAGCGTCGT"
        adapter = "AGATCGGAAGAGC"

        result = find_adapter_fuzzy(sequence, adapter)

        assert result is not None
        start, end = result
        assert start == 8
        assert sequence[start:end] == adapter

    def test_adapter_with_mismatch(self):
        """Test finding adapter with single mismatch."""
        sequence = "ATCGATCGAGATCGTAAGAGCGTCGT"  # G->T mismatch
        adapter = "AGATCGGAAGAGC"

        result = find_adapter_fuzzy(sequence, adapter, max_error_rate=0.15)

        assert result is not None
        # Should still find it with fuzzy matching

    def test_partial_adapter(self):
        """Test finding partial adapter at end."""
        sequence = "ATCGATCGAGATCGGAAG"
        adapter = "AGATCGGAAGAGC"  # Only first 10bp present

        result = find_adapter_fuzzy(sequence, adapter, min_overlap=6)

        assert result is not None
        # Should find partial match

    def test_no_adapter(self):
        """Test when adapter is not present."""
        sequence = "ATCGATCGATCGATCG"
        adapter = "AGATCGGAAGAGC"

        result = find_adapter_fuzzy(sequence, adapter)

        assert result is None

    def test_adapter_too_short(self):
        """Test with adapter shorter than minimum overlap."""
        sequence = "ATCGATCGATCGATCG"
        adapter = "AT"

        result = find_adapter_fuzzy(sequence, adapter, min_overlap=3)

        assert result is None

    def test_high_error_rate_rejected(self):
        """Test that high error rate alignments are rejected."""
        sequence = "ATCGATCGATCGATCG"
        adapter = "GGGGGGGGGGGG"  # No real match

        result = find_adapter_fuzzy(sequence, adapter, max_error_rate=0.1)

        # Should not match with low error tolerance
        assert result is None

    def test_adapter_at_start(self):
        """Test finding adapter at start of sequence."""
        adapter = "AGATCGGAAGAGC"
        sequence = adapter + "ATCGATCG"

        result = find_adapter_fuzzy(sequence, adapter)

        assert result is not None
        start, end = result
        assert start == 0
        assert sequence[start:end] == adapter


class TestCalculateMatchScore:
    """Test match score calculation."""

    def test_perfect_match(self):
        """Test perfect match score."""
        query = "ATCG"
        target = "ATCGNNNN"

        score = calculate_match_score(query, target, 0, 4)

        assert score == 1.0

    def test_partial_match(self):
        """Test partial match score."""
        query = "ATCG"
        target = "ATGGNNNN"  # 75% match (C->G)

        score = calculate_match_score(query, target, 0, 4)

        assert score == 0.75

    def test_no_match(self):
        """Test no match score."""
        query = "ATCG"
        target = "GGGGNNNN"

        score = calculate_match_score(query, target, 0, 4)

        assert score == 0.25  # Only G matches

    def test_invalid_positions(self):
        """Test with invalid positions."""
        query = "ATCG"
        target = "ATCG"

        # Start >= end
        score = calculate_match_score(query, target, 2, 2)
        assert score == 0.0

        # End > len(target)
        score = calculate_match_score(query, target, 0, 10)
        assert score == 0.0

    def test_case_insensitive(self):
        """Test case insensitive matching."""
        query = "ATCG"
        target = "atcgNNNN"

        score = calculate_match_score(query, target, 0, 4)

        assert score == 1.0
