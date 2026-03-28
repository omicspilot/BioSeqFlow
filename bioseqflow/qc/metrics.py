from __future__ import annotations

from bioseqflow.utils.io import FastqRecord

"""Custom quality metrics calculations."""

from typing import Any
from collections import defaultdict

import numpy as np

# Experiment-specific scoring profiles
SCORING_PROFILES = {
    "wgs": {  # Whole Genome Sequencing
        "base_quality": 0.35,
        "sequence_quality": 0.20,
        "length_distribution": 0.10,
        "gc_content": 0.15,
        "adapter_content": 0.15,
        "duplication": 0.05,  # Lower weight - some duplication OK
    },
    "rna_seq": {  # RNA Sequencing
        "base_quality": 0.30,
        "sequence_quality": 0.25,
        "length_distribution": 0.15,
        "gc_content": 0.25,
        "adapter_content": 0.05,
        "duplication": 0.00,  # Ignore - high duplication is normal
    },
    "amplicon": {  # Amplicon/Targeted Sequencing
        "base_quality": 0.40,
        "sequence_quality": 0.30,
        "length_distribution": 0.20,
        "gc_content": 0.10,
        "adapter_content": 0.00,
        "duplication": 0.00,  # Ignore - very high duplication expected
    },
    "chip_seq": {  # ChIP Sequencing
        "base_quality": 0.30,
        "sequence_quality": 0.20,
        "length_distribution": 0.10,
        "gc_content": 0.20,
        "adapter_content": 0.15,
        "duplication": 0.05,  # Some duplication expected at peaks
    },
}


class QualityMetrics:
    """Calculate custom quality metrics for sequencing data."""

    def __init__(self, experiment_type: str = "wgs") -> None:
        """
        Initialize quality metrics calculator.

        Args:
            experiment_type: Type of sequencing experiment
                ('wgs', 'rna_seq', 'amplicon', 'chip_seq')

        Raises:
            ValueError: If experiment_type is not recognized
        """
        if experiment_type not in SCORING_PROFILES:
            raise ValueError(
                f"Unknown experiment type: {experiment_type}. "
                f"Must be one of: {list(SCORING_PROFILES.keys())}"
            )
        self.experiment_type = experiment_type
        self.weights = SCORING_PROFILES[experiment_type]

    def calculate_composite_score(
        self, fastqc_data: dict[str, Any], _trimming_stats: dict[str, Any] | None = None
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
            _trimming_stats: Optional trimming statistics (unused currently)

        Returns:
            Composite quality score (0-100)
        """
        scores = {}

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

        BIOLOGICAL FIX: Context-aware scoring based on experiment type.
        High duplication is expected/normal for RNA-seq, amplicon, and ChIP-seq.

        Args:
            fastqc_data: FastQC data dictionary

        Returns:
            Score from 0-100
        """
        summary = fastqc_data.get("summary", {})
        module_status = summary.get("Sequence Duplication Levels", "FAIL")

        # Experiment-specific scoring
        if self.experiment_type in ["rna_seq", "amplicon"]:
            # Don't penalize expected high duplication
            return 100
        elif self.experiment_type == "chip_seq":
            # Moderate duplication acceptable (enriched regions)
            status_scores = {"PASS": 100, "WARN": 90, "FAIL": 70}
        else:  # wgs
            # Standard WGS - duplication is bad
            status_scores = {"PASS": 100, "WARN": 70, "FAIL": 30}

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

        CRITICAL FIX: Excludes N bases from calculation to avoid
        dilution of true GC content in low-quality regions.

        Args:
            sequence: DNA sequence string

        Returns:
            GC content as percentage (0-100)
        """
        sequence = sequence.upper()
        gc_count = sequence.count("G") + sequence.count("C")

        # Count only valid bases (exclude N and other ambiguity codes)
        valid_bases = sum(sequence.count(base) for base in "ACGT")

        if valid_bases == 0:
            return 0.0

        return (gc_count / valid_bases) * 100

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


    def calculate_per_base_quality(self, reads: list[FastqRecord], q_type: int = 33) -> dict[int, dict[str, float]]:
        """
        Calculate the statistics for each position of the different reads.
        Note:
            - mean_q: arithmetic mean of Q scores (what you have),
            - mean_error_prob: convert each Q → P, average P, convert back to Q

        Args:
            reads: a list of FastqRecord objects (varying lengths)
        Returns:
            a dict mapping each position index → a dict of stats (mean_q, mean_error_prob, median, q1, q3, min, max)
        """
        if not reads:
            return {}
        
        scores_by_position = defaultdict(list)
        for read in reads:
            # skip reads with no quality string
            if not read.quality:
                continue

            for pos in range(len(read.quality)):
                score = ord(read.quality[pos]) - q_type
                # append score in position i
                scores_by_position[pos].append(score)

        results = {}
        for pos in scores_by_position:
            scores = scores_by_position[pos]
            results[pos] = {
                "mean_q": float(np.mean(scores)),
                "median": float(np.median(scores)),
                "q1": float(np.percentile(scores, 25)),
                "q3": float(np.percentile(scores, 75)),
                "min": min(scores),
                "max": max(scores),
            }

        return results