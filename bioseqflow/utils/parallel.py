from __future__ import annotations

"""Parallel processing utilities for efficient batch operations."""

import multiprocessing as mp
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Any, Callable

from tqdm import tqdm


def get_optimal_workers(requested: int | None = None) -> int:
    """
    Get optimal number of worker processes.

    Args:
        requested: Requested number of workers (None = auto-detect)

    Returns:
        Number of workers to use
    """
    cpu_count = mp.cpu_count()

    if requested is None:
        # Use all cores minus one, minimum of 1
        return max(1, cpu_count - 1)

    # Ensure requested is within reasonable bounds
    return max(1, min(requested, cpu_count))


def process_samples_parallel(
    sample_list: Iterable[Any],
    function: Callable[[Any], Any],
    n_workers: int | None = None,
    use_threads: bool = False,
    show_progress: bool = True,
    description: str = "Processing",
) -> list[Any]:
    """
    Process multiple samples in parallel.

    Args:
        sample_list: Iterable of samples to process
        function: Function to apply to each sample
        n_workers: Number of parallel workers (None = auto-detect)
        use_threads: Use threads instead of processes (for I/O-bound tasks)
        show_progress: Show progress bar
        description: Progress bar description

    Returns:
        List of results in the same order as input

    Example:
        >>> def process_fastq(file_path):
        ...     return run_qc(file_path)
        >>> files = ["s1.fastq", "s2.fastq", "s3.fastq"]
        >>> results = process_samples_parallel(files, process_fastq, n_workers=4)
    """
    samples = list(sample_list)
    n_workers = get_optimal_workers(n_workers)

    executor_class = ThreadPoolExecutor if use_threads else ProcessPoolExecutor

    results = [None] * len(samples)

    with executor_class(max_workers=n_workers) as executor:
        # Submit all tasks
        future_to_idx = {
            executor.submit(function, sample): idx for idx, sample in enumerate(samples)
        }

        # Collect results with progress bar
        iterator = as_completed(future_to_idx)
        if show_progress:
            iterator = tqdm(
                iterator, total=len(samples), desc=description, unit="sample"
            )

        for future in iterator:
            idx = future_to_idx[future]
            try:
                results[idx] = future.result()
            except Exception as e:
                print(f"Error processing sample {idx}: {e}")
                results[idx] = None

    return results


def process_with_map(
    function: Callable[[Any], Any],
    items: Iterable[Any],
    n_workers: int | None = None,
    chunksize: int = 1,
) -> list[Any]:
    """
    Process items using multiprocessing.Pool.map.

    More memory efficient for large datasets than process_samples_parallel.

    Args:
        function: Function to apply
        items: Iterable of items to process
        n_workers: Number of workers
        chunksize: Number of items per chunk

    Returns:
        List of results
    """
    n_workers = get_optimal_workers(n_workers)

    with mp.Pool(processes=n_workers) as pool:
        results = pool.map(function, items, chunksize=chunksize)

    return results


def parallel_map(
    function: Callable[[Any], Any],
    items: Iterable[Any],
    n_workers: int | None = None,
) -> list[Any]:
    """
    Simple parallel map operation.

    Args:
        function: Function to apply
        items: Items to process
        n_workers: Number of workers

    Returns:
        List of results
    """
    return process_with_map(function, items, n_workers=n_workers)
