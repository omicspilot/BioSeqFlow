from __future__ import annotations

"""Custom quality metrics calculations."""

import statistics
from typing import Any

import numpy as np


class QualityMetrics:
    """Calculate custom quality metrics for sequencing data."""

    def __init__(self) -> None:
        """Initialize quality metrics calculator."""
        self.weights = {
            "base_quality": 0.3,
            "sequence_quality": 0.2,
            "length_distribution": 0.15,
            "gc_content": 0.15,
            "adapter_content": 0.1,
            "duplication": 0.1,
        }

    def calculate_composite_score(
        self, fastqc_data: dict[str, Any], trimming_stats: dict[str, Any] | None = None
    ) -> float:
        """
        Calculate composite quality score.

        Combines multiple quality metrics into a single score from 0-100:
        - >90: Excellent quality
        - 70-90: Good quality
        - 50-70: Needs review
        - <50: Requires attention

        Args:
            fastqc_data: Parsed FastQC data
            trimming_stats: Optional trimming statistics

        Returns:
            Composite quality score (0-100)
        """
        scores = {}

        # Extract basic statistics
        basic_stats = fastqc_data.get("modules", {}).get("Basic Statistics", {})

        # Base quality score (from per base sequence quality)
        scores["base_quality"] = self._score_base_quality(fastqc_data)

        # Sequence quality score
        scores["sequence_quality"] = self._score_sequence_quality(fastqc_data)

        # Length distribution score
        scores["length_distribution"] = self._score_length_distribution(fastqc_data)

        # GC content score
        scores["gc_content"] = self._score_gc_content(fastqc_data)

        # Adapter content score
        scores["adapter_content"] = self._score_adapter_content(fastqc_data)

        # Duplication score
        scores["duplication"] = self._score_duplication(fastqc_data)

        # Calculate weighted average
        composite_score = sum(
            scores[metric] * self.weights[metric] for metric in self.weights
        )

        return round(composite_score, 2)

    def _score_base_quality(self, fastqc_data: dict[str, Any]) -> float:
        """
        Score per-base sequence quality.

        Args:
            fastqc_data: FastQC data dictionary

        Returns:
            Score from 0-100
        """
        summary = fastqc_data.get("summary", {})
        module_status = summary.get("Per base sequence quality", "FAIL")

        status_scores = {"PASS": 100, "WARN": 70, "FAIL": 30}
        return status_scores.get(module_status, 50)

    def _score_sequence_quality(self, fastqc_data: dict[str, Any]) -> float:
        """
        Score per-sequence quality scores.

        Args:
            fastqc_data: FastQC data dictionary

        Returns:
            Score from 0-100
        """
        summary = fastqc_data.get("summary", {})
        module_status = summary.get("Per sequence quality scores", "FAIL")

        status_scores = {"PASS": 100, "WARN": 70, "FAIL": 30}
        return status_scores.get(module_status, 50)

    def _score_length_distribution(self, fastqc_data: dict[str, Any]) -> float:
        """
        Score sequence length distribution.

        Args:
            fastqc_data: FastQC data dictionary

        Returns:
            Score from 0-100
        """
        summary = fastqc_data.get("summary", {})
        module_status = summary.get("Sequence Length Distribution", "PASS")

        # Length distribution warnings are often acceptable
        status_scores = {"PASS": 100, "WARN": 85, "FAIL": 60}
        return status_scores.get(module_status, 80)

    def _score_gc_content(self, fastqc_data: dict[str, Any]) -> float:
        """
        Score GC content distribution.

        Args:
            fastqc_data: FastQC data dictionary

        Returns:
            Score from 0-100
        """
        summary = fastqc_data.get("summary", {})
        module_status = summary.get("Per sequence GC content", "FAIL")

        status_scores = {"PASS": 100, "WARN": 75, "FAIL": 40}
        return status_scores.get(module_status, 50)

    def _score_adapter_content(self, fastqc_data: dict[str, Any]) -> float:
        """
        Score adapter content.

        Args:
            fastqc_data: FastQC data dictionary

        Returns:
            Score from 0-100
        """
        summary = fastqc_data.get("summary", {})
        module_status = summary.get("Adapter Content", "PASS")

        status_scores = {"PASS": 100, "WARN": 60, "FAIL": 20}
        return status_scores.get(module_status, 50)

    def _score_duplication(self, fastqc_data: dict[str, Any]) -> float:
        """
        Score sequence duplication levels.

        Args:
            fastqc_data: FastQC data dictionary

        Returns:
            Score from 0-100
        """
        summary = fastqc_data.get("summary", {})
        module_status = summary.get("Sequence Duplication Levels", "FAIL")

        # High duplication can be normal for RNA-seq or amplicon sequencing
        status_scores = {"PASS": 100, "WARN": 80, "FAIL": 50}
        return status_scores.get(module_status, 60)

    def calculate_quality_distribution(
        self, quality_scores: list[int]
    ) -> dict[str, float]:
        """
        Calculate quality score distribution statistics.

        Args:
            quality_scores: List of Phred quality scores

        Returns:
            Dictionary with mean, median, std, min, max, and quartiles
        """
        if not quality_scores:
            return {}

        quality_array = np.array(quality_scores)

        return {
            "mean": float(np.mean(quality_array)),
            "median": float(np.median(quality_array)),
            "std": float(np.std(quality_array)),
            "min": float(np.min(quality_array)),
            "max": float(np.max(quality_array)),
            "q25": float(np.percentile(quality_array, 25)),
            "q75": float(np.percentile(quality_array, 75)),
        }

    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate GC content percentage.

        Args:
            sequence: DNA sequence string

        Returns:
            GC content as percentage (0-100)
        """
        sequence = sequence.upper()
        gc_count = sequence.count("G") + sequence.count("C")
        total = len(sequence)

        if total == 0:
            return 0.0

        return (gc_count / total) * 100

    def calculate_n_content(self, sequence: str) -> float:
        """
        Calculate N content percentage.

        Args:
            sequence: DNA sequence string

        Returns:
            N content as percentage (0-100)
        """
        sequence = sequence.upper()
        n_count = sequence.count("N")
        total = len(sequence)

        if total == 0:
            return 0.0

        return (n_count / total) * 100
