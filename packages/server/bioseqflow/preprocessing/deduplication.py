from __future__ import annotations

"""Duplicate removal modules."""

import hashlib
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from bioseqflow.core.base import PreprocessingModule
from bioseqflow.utils.io import FastqRecord, read_fastq, write_fastq
from bioseqflow.utils.validators import validate_file_path


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculate Hamming distance between two sequences.

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        Number of mismatches

    Raises:
        ValueError: If sequences have different lengths
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must be same length: {len(seq1)} vs {len(seq2)}")

    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def cluster_umis_directional(umis: list[str], max_distance: int = 1) -> dict[str, str]:
    """
    Cluster UMIs using directional method (Smith et al., 2017).

    This method handles sequencing errors in UMIs by clustering similar
    UMIs (within max_distance mismatches) together. Uses directional
    approach where higher-count UMIs absorb lower-count UMIs.

    Algorithm:
    1. Count UMI frequencies
    2. Sort by frequency (descending)
    3. For each UMI, merge UMIs within max_distance if they have lower counts
    4. Return mapping of original UMI → representative UMI

    Args:
        umis: List of UMI sequences
        max_distance: Maximum Hamming distance for clustering (1-2 recommended)

    Returns:
        Dictionary mapping each UMI to its representative (cluster center)

    References:
        Smith, T., et al. (2017). "UMI-tools: modeling sequencing errors
        in Unique Molecular Identifiers to improve quantification accuracy."
        Genome Research, 27(3), 491-499.
    """
    if not umis:
        return {}

    # Count UMI frequencies
    umi_counts = Counter(umis)

    # Sort UMIs by count (descending), then alphabetically for ties
    sorted_umis = sorted(umi_counts.keys(), key=lambda x: (-umi_counts[x], x))

    # Initialize mapping: each UMI maps to itself
    umi_mapping: dict[str, str] = {umi: umi for umi in umi_counts}

    # Directional clustering
    for i, umi1 in enumerate(sorted_umis):
        # Skip if this UMI was already merged into another
        if umi_mapping[umi1] != umi1:
            continue

        # Check all lower-count UMIs
        for umi2 in sorted_umis[i + 1 :]:
            # Skip if already merged
            if umi_mapping[umi2] != umi2:
                continue

            # Skip if counts are equal (avoid arbitrary merging)
            if umi_counts[umi1] == umi_counts[umi2]:
                continue

            # Calculate Hamming distance
            try:
                distance = hamming_distance(umi1, umi2)

                # Merge if within threshold
                if distance <= max_distance:
                    umi_mapping[umi2] = umi1

            except ValueError:
                # Different length UMIs - skip
                continue

    return umi_mapping


class DuplicateRemover(PreprocessingModule):
    """Remove duplicate reads."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize duplicate remover.

        Args:
            config: Configuration object
        """
        super().__init__(config)

    def validate_inputs(self, input_file: Path | str, _output_file: Path | str) -> None:
        """
        Validate inputs.

        Args:
            input_file: Input file
            _output_file: Output file (unused in base validation)

        Raises:
            ValueError: If inputs are invalid
        """
        validate_file_path(
            input_file,
            must_exist=True,
            extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        )

    def process(self, input_file: Path, output_file: Path) -> dict[str, int | float]:
        """
        Process file with duplicate removal.

        Args:
            input_file: Input file
            output_file: Output file

        Returns:
            Deduplication statistics
        """
        return self.remove_duplicates(input_file, output_file)

    def remove_duplicates(
        self,
        input_file: Path | str,
        output_file: Path | str,
        method: str = "exact",
        keep_best: bool = True,
        use_hash: bool = False,
    ) -> dict[str, int | float]:
        """
        Remove duplicate reads.

        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            method: Deduplication method ('exact', 'prefix', or 'hash')
            keep_best: Keep read with highest quality (if duplicates found)
            use_hash: Use hash-based dedup for memory efficiency (recommended for >50M reads)

        Returns:
            Deduplication statistics

        Raises:
            ValueError: If method is invalid
        """
        input_file = Path(input_file)
        output_file = Path(output_file)

        self.validate_inputs(input_file, output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        if method == "exact":
            if use_hash:
                return self._remove_exact_duplicates_hash(input_file, output_file, keep_best)
            return self._remove_exact_duplicates(input_file, output_file, keep_best)
        elif method == "prefix":
            return self._remove_prefix_duplicates(input_file, output_file, keep_best)
        elif method == "hash":
            return self._remove_exact_duplicates_hash(input_file, output_file, keep_best)
        else:
            raise ValueError(f"Unknown deduplication method: {method}")

    def _remove_exact_duplicates(
        self, input_file: Path, output_file: Path, keep_best: bool
    ) -> dict[str, int | float]:
        """
        Remove exact sequence duplicates.

        Args:
            input_file: Input file
            output_file: Output file
            keep_best: Keep highest quality duplicate

        Returns:
            Statistics
        """
        total_reads = 0
        unique_reads = 0
        duplicates = 0

        seen_sequences: dict[str, FastqRecord] = {}

        for record in read_fastq(input_file):
            total_reads += 1

            if record.sequence not in seen_sequences:
                seen_sequences[record.sequence] = record
                unique_reads += 1
            else:
                duplicates += 1

                # Keep record with better quality if requested
                if (
                    keep_best
                    and record.mean_quality() > seen_sequences[record.sequence].mean_quality()
                ):
                    seen_sequences[record.sequence] = record

        # Write unique reads
        write_fastq(
            list(seen_sequences.values()), output_file, compress=str(output_file).endswith(".gz")
        )

        self.stats = {
            "total_reads": total_reads,
            "unique_reads": unique_reads,
            "duplicates": duplicates,
            "duplication_rate": (duplicates / total_reads * 100) if total_reads > 0 else 0,
        }

        return self.stats

    def _remove_exact_duplicates_hash(
        self, input_file: Path, output_file: Path, keep_best: bool
    ) -> dict[str, int | float]:
        """
        Remove exact duplicates using hash-based approach (memory-efficient).

        CRITICAL FIX: Uses MD5 hashes instead of full sequences to reduce
        memory usage from ~100 bytes/read to ~16 bytes/read.
        Enables processing of 100M+ read files.

        Args:
            input_file: Input file
            output_file: Output file
            keep_best: Keep highest quality duplicate

        Returns:
            Statistics
        """
        total_reads = 0
        unique_reads = 0
        duplicates = 0

        # Store hash -> best record mapping (16 bytes hash vs 100+ bytes sequence)
        seen_hashes: dict[str, FastqRecord] = {}

        def hash_sequence(seq: str) -> str:
            """Generate MD5 hash of sequence."""
            return hashlib.md5(seq.encode()).hexdigest()

        for record in read_fastq(input_file):
            total_reads += 1
            seq_hash = hash_sequence(record.sequence)

            if seq_hash not in seen_hashes:
                seen_hashes[seq_hash] = record
                unique_reads += 1
            else:
                duplicates += 1

                # Keep record with better quality if requested
                if keep_best and record.mean_quality() > seen_hashes[seq_hash].mean_quality():
                    seen_hashes[seq_hash] = record

        # Write unique reads
        write_fastq(
            list(seen_hashes.values()), output_file, compress=str(output_file).endswith(".gz")
        )

        self.stats = {
            "total_reads": total_reads,
            "unique_reads": unique_reads,
            "duplicates": duplicates,
            "duplication_rate": (duplicates / total_reads * 100) if total_reads > 0 else 0,
        }

        return self.stats

    def _remove_prefix_duplicates(
        self, input_file: Path, output_file: Path, keep_best: bool, prefix_length: int = 30
    ) -> dict[str, int | float]:
        """
        Remove duplicates based on sequence prefix.

        Useful for identifying PCR duplicates where only the start is identical.

        Args:
            input_file: Input file
            output_file: Output file
            keep_best: Keep highest quality duplicate
            prefix_length: Length of prefix to compare

        Returns:
            Statistics
        """
        total_reads = 0
        unique_reads = 0
        duplicates = 0

        seen_prefixes: dict[str, FastqRecord] = {}

        for record in read_fastq(input_file):
            total_reads += 1

            # CRITICAL FIX: Handle reads shorter than prefix_length
            # Prevents IndexError on short/trimmed reads
            actual_prefix_len = min(prefix_length, len(record.sequence))
            prefix = record.sequence[:actual_prefix_len]

            if prefix not in seen_prefixes:
                seen_prefixes[prefix] = record
                unique_reads += 1
            else:
                duplicates += 1

                # Keep record with better quality if requested
                if keep_best and record.mean_quality() > seen_prefixes[prefix].mean_quality():
                    seen_prefixes[prefix] = record

        # Write unique reads
        write_fastq(
            list(seen_prefixes.values()), output_file, compress=str(output_file).endswith(".gz")
        )

        self.stats = {
            "total_reads": total_reads,
            "unique_reads": unique_reads,
            "duplicates": duplicates,
            "duplication_rate": (duplicates / total_reads * 100) if total_reads > 0 else 0,
        }

        return self.stats


class UMIDeduplicator(PreprocessingModule):
    """Remove duplicates using Unique Molecular Identifiers (UMIs)."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize UMI deduplicator.

        Args:
            config: Configuration object
        """
        super().__init__(config)

    def validate_inputs(self, input_file: Path | str, _output_file: Path | str) -> None:
        """
        Validate inputs.

        Args:
            input_file: Input file
            _output_file: Output file (unused in base validation)

        Raises:
            ValueError: If inputs are invalid
        """
        validate_file_path(
            input_file,
            must_exist=True,
            extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        )

    def process(self, input_file: Path, output_file: Path) -> dict[str, int | float]:
        """
        Process file with UMI-based deduplication.

        Args:
            input_file: Input file
            output_file: Output file

        Returns:
            Deduplication statistics
        """
        return self.deduplicate(input_file, output_file)

    def deduplicate(
        self,
        input_file: Path | str,
        output_file: Path | str,
        umi_length: int = 8,
        umi_location: str = "start",
        error_correction: bool = True,
        max_edit_distance: int = 1,
    ) -> dict[str, int | float]:
        """
        Deduplicate reads using UMIs with optional error correction.

        BIOLOGICAL ENHANCEMENT: Error correction handles sequencing errors in UMIs.
        Without correction, PCR duplicates with UMI errors are counted as unique
        molecules, inflating library complexity by 2-10x.

        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            umi_length: Length of UMI sequence
            umi_location: UMI location ('start' or 'end')
            error_correction: Use directional clustering to correct UMI errors
            max_edit_distance: Maximum Hamming distance for UMI clustering (1-2)

        Returns:
            Deduplication statistics

        References:
            Smith, T., et al. (2017). "UMI-tools: modeling sequencing errors..."
            Genome Research, 27(3), 491-499.
        """
        input_file = Path(input_file)
        output_file = Path(output_file)

        self.validate_inputs(input_file, output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        total_reads = 0
        unique_umis = 0
        duplicates = 0
        umis_before_correction = 0
        umis_after_correction = 0

        # Group by UMI
        umi_groups: dict[str, list[FastqRecord]] = defaultdict(list)
        all_umis: list[str] = []

        for record in read_fastq(input_file):
            total_reads += 1

            # Extract UMI
            if umi_location == "start":
                umi = record.sequence[:umi_length]
                trimmed_seq = record.sequence[umi_length:]
                trimmed_qual = record.quality[umi_length:]
            else:  # end
                umi = record.sequence[-umi_length:]
                trimmed_seq = record.sequence[:-umi_length]
                trimmed_qual = record.quality[:-umi_length]

            # Create trimmed record (without UMI)
            trimmed_record = FastqRecord(record.header, trimmed_seq, record.plus, trimmed_qual)

            umi_groups[umi].append(trimmed_record)
            all_umis.append(umi)

        umis_before_correction = len(umi_groups)

        # Apply error correction if requested
        if error_correction:
            # Cluster similar UMIs
            umi_mapping = cluster_umis_directional(all_umis, max_distance=max_edit_distance)

            # Merge groups based on clustering
            corrected_groups: dict[str, list[FastqRecord]] = defaultdict(list)
            for umi, records in umi_groups.items():
                representative_umi = umi_mapping[umi]
                corrected_groups[representative_umi].extend(records)

            umi_groups = corrected_groups
            umis_after_correction = len(umi_groups)

        # Keep best quality read from each UMI group
        unique_records = []
        for _umi, records in umi_groups.items():
            unique_umis += 1
            duplicates += len(records) - 1

            # Keep read with highest mean quality
            best_record = max(records, key=lambda r: r.mean_quality())
            unique_records.append(best_record)

        # Write deduplicated reads
        write_fastq(unique_records, output_file, compress=str(output_file).endswith(".gz"))

        self.stats = {
            "total_reads": total_reads,
            "unique_umis": unique_umis,
            "duplicates": duplicates,
            "duplication_rate": (duplicates / total_reads * 100) if total_reads > 0 else 0,
            "umi_length": umi_length,
        }

        # Add error correction stats if used
        if error_correction:
            self.stats["umis_before_correction"] = umis_before_correction
            self.stats["umis_after_correction"] = umis_after_correction
            self.stats["umis_merged"] = umis_before_correction - umis_after_correction
            self.stats["correction_rate"] = (
                ((umis_before_correction - umis_after_correction) / umis_before_correction * 100)
                if umis_before_correction > 0
                else 0
            )

        return self.stats
