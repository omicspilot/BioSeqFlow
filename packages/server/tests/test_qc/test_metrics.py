"""Tests for quality metrics."""

import pytest

from bioseqflow_server.toolkit.qc.metrics import QualityMetrics
from bioseqflow_server.toolkit.utils.io import FastqRecord


class TestQualityMetrics:
    """Tests for QualityMetrics class."""

    @pytest.fixture
    def metrics(self) -> QualityMetrics:
        """Create QualityMetrics instance."""
        return QualityMetrics()

    @pytest.fixture
    def sample_fastqc_data(self) -> dict:  # type: ignore[type-arg]
        """Sample FastQC data for testing."""
        return {
            "summary": {
                "Per base sequence quality": "PASS",
                "Per sequence quality scores": "PASS",
                "Sequence Length Distribution": "PASS",
                "Per sequence GC content": "PASS",
                "Adapter Content": "PASS",
                "Sequence Duplication Levels": "WARN",
            },
            "modules": {
                "Basic Statistics": {
                    "Total Sequences": 1000000,
                    "%GC": 48,
                }
            },
        }

    def test_calculate_composite_score(
        self,
        metrics: QualityMetrics,
        sample_fastqc_data: dict,  # type: ignore[type-arg]
    ) -> None:
        """Test composite score calculation."""
        score = metrics.calculate_composite_score(sample_fastqc_data)

        assert isinstance(score, float)
        assert 0 <= score <= 100

    def test_score_base_quality_pass(
        self,
        metrics: QualityMetrics,
        sample_fastqc_data: dict,  # type: ignore[type-arg]
    ) -> None:
        """Test base quality scoring with PASS."""
        score = metrics._score_base_quality(sample_fastqc_data)

        assert score == 100

    def test_score_base_quality_warn(
        self,
        metrics: QualityMetrics,
        sample_fastqc_data: dict,  # type: ignore[type-arg]
    ) -> None:
        """Test base quality scoring with WARN."""
        sample_fastqc_data["summary"]["Per base sequence quality"] = "WARN"
        score = metrics._score_base_quality(sample_fastqc_data)

        assert score == 70

    def test_score_base_quality_fail(
        self,
        metrics: QualityMetrics,
        sample_fastqc_data: dict,  # type: ignore[type-arg]
    ) -> None:
        """Test base quality scoring with FAIL."""
        sample_fastqc_data["summary"]["Per base sequence quality"] = "FAIL"
        score = metrics._score_base_quality(sample_fastqc_data)

        assert score == 30

    def test_calculate_quality_distribution(self, metrics: QualityMetrics) -> None:
        """Test quality distribution calculation."""
        quality_scores = [20, 25, 30, 35, 40, 30, 25, 30]
        dist = metrics.calculate_quality_distribution(quality_scores)

        assert "mean" in dist
        assert "median" in dist
        assert "std" in dist
        assert "min" in dist
        assert "max" in dist

        assert dist["min"] == 20
        assert dist["max"] == 40

    def test_calculate_quality_distribution_empty(self, metrics: QualityMetrics) -> None:
        """Test quality distribution with empty list."""
        dist = metrics.calculate_quality_distribution([])

        assert dist == {}

    def test_calculate_gc_content(self, metrics: QualityMetrics) -> None:
        """Test GC content calculation."""
        # 50% GC (4 G/C out of 8 bases)
        sequence = "ATCGATCG"
        gc_content = metrics.calculate_gc_content(sequence)

        assert gc_content == 50.0

    def test_calculate_gc_content_all_gc(self, metrics: QualityMetrics) -> None:
        """Test GC content with all G/C."""
        sequence = "GCGCGCGC"
        gc_content = metrics.calculate_gc_content(sequence)

        assert gc_content == 100.0

    def test_calculate_gc_content_no_gc(self, metrics: QualityMetrics) -> None:
        """Test GC content with no G/C."""
        sequence = "ATATATAT"
        gc_content = metrics.calculate_gc_content(sequence)

        assert gc_content == 0.0

    def test_calculate_gc_content_empty(self, metrics: QualityMetrics) -> None:
        """Test GC content with empty sequence."""
        gc_content = metrics.calculate_gc_content("")

        assert gc_content == 0.0

    def test_calculate_n_content(self, metrics: QualityMetrics) -> None:
        """Test N content calculation."""
        sequence = "ATCGNNATCG"
        n_content = metrics.calculate_n_content(sequence)

        assert n_content == 20.0  # 2 N's out of 10 bases

    def test_calculate_n_content_no_n(self, metrics: QualityMetrics) -> None:
        """Test N content with no N's."""
        sequence = "ATCGATCG"
        n_content = metrics.calculate_n_content(sequence)

        assert n_content == 0.0

    def test_per_base_quality_known_values(self, metrics: QualityMetrics) -> None:
        """
        Stats at each position are correct for known quality strings.
        """
        # 'I' = ASCII 73, Q = 73-33 = 40
        # '5' = ASCII 53, Q = 53-33 = 20
        # Two reads, same length, deterministic values
        reads = [
            FastqRecord("@read1", "ACGT", "+", "IIII"),  # all Q40
            FastqRecord("@read2", "ACGT", "+", "5555"),  # all Q20
        ]
        result = metrics.calculate_per_base_quality(reads)

        # Position 0: scores are [40, 20] → mean=30, min=20, max=40
        assert result[0]["mean_q"] == 30.0
        assert result[0]["min"] == 20
        assert result[0]["max"] == 40

    def test_per_base_quality_empty_reads(self, metrics: QualityMetrics) -> None:
        """Empty input returns empty dict."""
        result = metrics.calculate_per_base_quality([])
        assert result == {}
