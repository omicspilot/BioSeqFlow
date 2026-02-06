from __future__ import annotations

"""Tests for parallel processing utilities."""

import multiprocessing as mp

import pytest

from bioseqflow.utils.parallel import (
    get_optimal_workers,
    parallel_map,
    process_samples_parallel,
    process_with_map,
)


def simple_function(x: int) -> int:
    """Simple function for testing - doubles the input."""
    return x * 2


def failing_function(x: int) -> int:
    """Function that raises an error for testing error handling."""
    if x < 0:
        raise ValueError("Negative value not allowed")
    return x * 2


class TestGetOptimalWorkers:
    """Test optimal worker calculation."""

    def test_get_optimal_workers_auto(self):
        """Test automatic worker detection."""
        workers = get_optimal_workers(None)
        cpu_count = mp.cpu_count()

        # Should be all cores minus one, minimum 1
        expected = max(1, cpu_count - 1)
        assert workers == expected

    def test_get_optimal_workers_requested(self):
        """Test with requested number of workers."""
        workers = get_optimal_workers(4)
        assert workers == min(4, mp.cpu_count())

    def test_get_optimal_workers_exceeds_cpu(self):
        """Test when requested workers exceed CPU count."""
        cpu_count = mp.cpu_count()
        workers = get_optimal_workers(cpu_count + 10)

        # Should cap at CPU count
        assert workers == cpu_count

    def test_get_optimal_workers_zero(self):
        """Test with zero workers requested."""
        workers = get_optimal_workers(0)

        # Should use at least 1 worker
        assert workers == 1

    def test_get_optimal_workers_negative(self):
        """Test with negative workers requested."""
        workers = get_optimal_workers(-5)

        # Should use at least 1 worker
        assert workers == 1


class TestProcessSamplesParallel:
    """Test parallel sample processing."""

    def test_process_samples_parallel_basic(self):
        """Test basic parallel processing."""
        samples = [1, 2, 3, 4, 5]
        results = process_samples_parallel(
            samples, simple_function, n_workers=2, show_progress=False
        )

        assert results == [2, 4, 6, 8, 10]

    def test_process_samples_parallel_with_threads(self):
        """Test parallel processing with threads."""
        samples = [1, 2, 3]
        results = process_samples_parallel(
            samples, simple_function, n_workers=2, use_threads=True, show_progress=False
        )

        assert results == [2, 4, 6]

    def test_process_samples_parallel_auto_workers(self):
        """Test with automatic worker detection."""
        samples = [1, 2, 3]
        results = process_samples_parallel(
            samples, simple_function, n_workers=None, show_progress=False
        )

        assert results == [2, 4, 6]

    def test_process_samples_parallel_empty_list(self):
        """Test with empty sample list."""
        samples: list[int] = []
        results = process_samples_parallel(
            samples, simple_function, show_progress=False
        )

        assert results == []

    def test_process_samples_parallel_single_item(self):
        """Test with single item."""
        samples = [5]
        results = process_samples_parallel(
            samples, simple_function, show_progress=False
        )

        assert results == [10]

    def test_process_samples_parallel_with_error(self):
        """Test error handling in parallel processing."""
        samples = [1, -1, 3]  # -1 will cause error
        results = process_samples_parallel(
            samples, failing_function, n_workers=2, show_progress=False
        )

        # Should have None for failed item
        assert results[0] == 2
        assert results[1] is None  # Failed
        assert results[2] == 6

    def test_process_samples_parallel_order_preserved(self):
        """Test that output order matches input order."""
        samples = list(range(10))
        results = process_samples_parallel(
            samples, simple_function, n_workers=3, show_progress=False
        )

        expected = [x * 2 for x in samples]
        assert results == expected

    def test_process_samples_parallel_with_progress(self):
        """Test parallel processing with progress bar."""
        samples = [1, 2, 3]
        # Just test that it doesn't crash with progress enabled
        results = process_samples_parallel(
            samples,
            simple_function,
            n_workers=2,
            show_progress=True,
            description="Testing"
        )

        assert results == [2, 4, 6]


class TestProcessWithMap:
    """Test process_with_map function."""

    def test_process_with_map_basic(self):
        """Test basic map processing."""
        items = [1, 2, 3, 4, 5]
        results = process_with_map(simple_function, items, n_workers=2)

        assert results == [2, 4, 6, 8, 10]

    def test_process_with_map_auto_workers(self):
        """Test with automatic worker detection."""
        items = [1, 2, 3]
        results = process_with_map(simple_function, items, n_workers=None)

        assert results == [2, 4, 6]

    def test_process_with_map_chunksize(self):
        """Test with custom chunksize."""
        items = list(range(10))
        results = process_with_map(
            simple_function, items, n_workers=2, chunksize=2
        )

        expected = [x * 2 for x in range(10)]
        assert results == expected

    def test_process_with_map_empty_list(self):
        """Test with empty list."""
        items: list[int] = []
        results = process_with_map(simple_function, items)

        assert results == []

    def test_process_with_map_single_item(self):
        """Test with single item."""
        items = [7]
        results = process_with_map(simple_function, items)

        assert results == [14]


class TestParallelMap:
    """Test parallel_map function."""

    def test_parallel_map_basic(self):
        """Test basic parallel map."""
        items = [1, 2, 3, 4]
        results = parallel_map(simple_function, items, n_workers=2)

        assert results == [2, 4, 6, 8]

    def test_parallel_map_auto_workers(self):
        """Test with automatic worker detection."""
        items = [1, 2, 3]
        results = parallel_map(simple_function, items, n_workers=None)

        assert results == [2, 4, 6]

    def test_parallel_map_order_preserved(self):
        """Test that order is preserved."""
        items = list(range(20))
        results = parallel_map(simple_function, items, n_workers=4)

        expected = [x * 2 for x in range(20)]
        assert results == expected
