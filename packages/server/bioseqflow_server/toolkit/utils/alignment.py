from __future__ import annotations

"""Sequence alignment utilities for adapter matching."""

from typing import NamedTuple


class AlignmentResult(NamedTuple):
    """Result from sequence alignment."""

    score: int
    start: int
    end: int
    aligned_query: str
    aligned_target: str


def smith_waterman(
    query: str,
    target: str,
    match_score: int = 2,
    mismatch_penalty: int = -1,
    gap_penalty: int = -1,
) -> AlignmentResult:
    """
    Perform Smith-Waterman local alignment.

    This algorithm finds the best local alignment between query and target
    sequences, allowing mismatches, insertions, and deletions. Ideal for
    fuzzy adapter matching where adapters may have sequencing errors or
    be partially present.

    Args:
        query: Query sequence (e.g., adapter)
        target: Target sequence (e.g., read)
        match_score: Score for matching bases
        mismatch_penalty: Penalty for mismatched bases
        gap_penalty: Penalty for gaps (insertions/deletions)

    Returns:
        AlignmentResult with score and alignment positions

    References:
        Smith, T.F. and Waterman, M.S. (1981). "Identification of common
        molecular subsequences." Journal of Molecular Biology, 147(1), 195-197.
    """
    m, n = len(query), len(target)

    # Initialize scoring matrix (m+1 x n+1)
    matrix = [[0] * (n + 1) for _ in range(m + 1)]

    # Track maximum score position
    max_score = 0
    max_pos = (0, 0)

    # Fill scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Score for match/mismatch
            if query[i - 1].upper() == target[j - 1].upper():
                diagonal = matrix[i - 1][j - 1] + match_score
            else:
                diagonal = matrix[i - 1][j - 1] + mismatch_penalty

            # Score for gaps
            up = matrix[i - 1][j] + gap_penalty
            left = matrix[i][j - 1] + gap_penalty

            # Take maximum (or 0 for local alignment)
            matrix[i][j] = max(0, diagonal, up, left)

            # Track maximum score
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    # Traceback to find alignment
    i, j = max_pos
    aligned_query: list[str] = []
    aligned_target: list[str] = []

    # Store end position in target
    end_pos = j

    # Traceback until we hit 0
    while i > 0 and j > 0 and matrix[i][j] > 0:
        current = matrix[i][j]
        diagonal = matrix[i - 1][j - 1]
        up = matrix[i - 1][j]
        left = matrix[i][j - 1]

        if query[i - 1].upper() == target[j - 1].upper():
            diagonal_score = diagonal + match_score
        else:
            diagonal_score = diagonal + mismatch_penalty

        if current == diagonal_score:
            # Match/mismatch
            aligned_query.insert(0, query[i - 1])
            aligned_target.insert(0, target[j - 1])
            i -= 1
            j -= 1
        elif current == up + gap_penalty:
            # Gap in target
            aligned_query.insert(0, query[i - 1])
            aligned_target.insert(0, "-")
            i -= 1
        else:
            # Gap in query
            aligned_query.insert(0, "-")
            aligned_target.insert(0, target[j - 1])
            j -= 1

    # Start position in target
    start_pos = j

    return AlignmentResult(
        score=max_score,
        start=start_pos,
        end=end_pos,
        aligned_query="".join(aligned_query),
        aligned_target="".join(aligned_target),
    )


def find_adapter_fuzzy(
    sequence: str,
    adapter: str,
    min_overlap: int = 3,
    max_error_rate: float = 0.15,
) -> tuple[int, int] | None:
    """
    Find adapter in sequence using fuzzy matching.

    Uses Smith-Waterman alignment to locate adapter allowing for
    mismatches and partial matches.

    Args:
        sequence: Read sequence
        adapter: Adapter sequence to find
        min_overlap: Minimum overlap length to consider
        max_error_rate: Maximum error rate (mismatches/length)

    Returns:
        Tuple of (start, end) positions if found, None otherwise

    Examples:
        >>> find_adapter_fuzzy("ATCGATCGAGATCGGAAGAG", "AGATCGGAAGAGC")
        (8, 20)  # Found with partial match

        >>> find_adapter_fuzzy("ATCGATCG", "AGATCGGAAGAGC")
        None  # No significant match
    """
    if len(adapter) < min_overlap:
        return None

    # Perform Smith-Waterman alignment
    alignment = smith_waterman(adapter, sequence)

    # Check if alignment meets minimum criteria
    aligned_length = len(alignment.aligned_query.replace("-", ""))

    if aligned_length < min_overlap:
        return None

    # Calculate error rate (mismatches + gaps)
    mismatches = sum(1 for q, t in zip(alignment.aligned_query, alignment.aligned_target) if q != t)
    error_rate = mismatches / aligned_length if aligned_length > 0 else 1.0

    if error_rate > max_error_rate:
        return None

    return (alignment.start, alignment.end)


def calculate_match_score(query: str, target: str, start: int, end: int) -> float:
    """
    Calculate alignment quality score.

    Args:
        query: Query sequence (adapter)
        target: Target sequence (read)
        start: Start position in target
        end: End position in target

    Returns:
        Match score (0.0 to 1.0)
    """
    if start >= end or end > len(target):
        return 0.0

    aligned_region = target[start:end]
    query_prefix = query[: len(aligned_region)]

    if len(query_prefix) == 0:
        return 0.0

    matches = sum(1 for q, t in zip(query_prefix.upper(), aligned_region.upper()) if q == t)

    return matches / len(query_prefix)
