#!/usr/bin/env python3
"""Example: Parallel processing for batch QC analysis.

This example demonstrates:
1. Automatic worker detection
2. Process-based parallel processing
3. Thread-based parallel processing
4. Different parallel processing strategies
5. Error handling in parallel workflows
"""

import time
from pathlib import Path

from bioseqflow.preprocessing.filtering import QualityFilter
from bioseqflow.utils.io import FastqRecord, write_fastq
from bioseqflow.utils.parallel import (
    get_optimal_workers,
    parallel_map,
    process_samples_parallel,
    process_with_map,
)


def create_sample_files(output_dir: Path, n_samples: int = 5) -> list[Path]:
    """Create sample FASTQ files for demonstration."""
    output_dir.mkdir(exist_ok=True)

    sample_files = []
    for i in range(1, n_samples + 1):
        sample_file = output_dir / f"sample_{i}.fastq"

        # Create records with varying quality
        records = []
        for j in range(100):
            seq = "ATCGATCG" * 10  # 80bp
            # Vary quality: samples get progressively better quality
            qual_char = chr(33 + 20 + i * 5)  # Q20, Q25, Q30, etc.
            qual = qual_char * len(seq)
            records.append(
                FastqRecord(
                    f"@read_{j}\n",
                    seq + "\n",
                    "+\n",
                    qual + "\n",
                )
            )

        write_fastq(records, sample_file)
        sample_files.append(sample_file)

    return sample_files


def process_single_sample(sample_file: Path) -> dict[str, int | float]:
    """Process a single FASTQ file (simulates QC analysis)."""
    print(f"  Processing {sample_file.name}...")

    # Simulate some processing time
    time.sleep(0.5)

    # Run quality filtering
    output_file = sample_file.parent / f"{sample_file.stem}_filtered.fastq"
    filter_obj = QualityFilter()

    stats = filter_obj.filter(
        sample_file,
        output_file,
        min_quality=25,
        min_length=50
    )

    return {
        "sample": sample_file.name,
        "total_reads": stats["total_reads"],
        "passed_reads": stats["passed_reads"],
        "pass_rate": (stats["passed_reads"] / stats["total_reads"] * 100)
        if stats["total_reads"] > 0 else 0
    }


def example_worker_detection():
    """Demonstrate automatic worker detection."""
    print("=" * 60)
    print("Worker Detection")
    print("=" * 60)

    # Get optimal workers (auto-detect)
    optimal = get_optimal_workers(None)
    print(f"\nOptimal workers (auto): {optimal}")

    # Request specific number
    requested = get_optimal_workers(4)
    print(f"Requested 4 workers: {requested}")

    # Request more than available
    high = get_optimal_workers(1000)
    print(f"Requested 1000 workers: {high} (capped at CPU count)")

    # Request zero/negative (gets capped to 1)
    low = get_optimal_workers(0)
    print(f"Requested 0 workers: {low} (minimum is 1)")


def example_process_samples_parallel(sample_files: list[Path]):
    """Demonstrate process_samples_parallel with progress bar."""
    print("\n" + "=" * 60)
    print("Method 1: process_samples_parallel (with progress bar)")
    print("=" * 60)

    print("\nProcessing samples in parallel...")
    start_time = time.time()

    results = process_samples_parallel(
        sample_files,
        process_single_sample,
        n_workers=3,
        use_threads=False,
        show_progress=True,
        description="QC Analysis"
    )

    elapsed = time.time() - start_time

    print(f"\nCompleted in {elapsed:.2f} seconds")
    print("\nResults:")
    for result in results:
        if result:
            print(f"  {result['sample']}: "
                  f"{result['passed_reads']}/{result['total_reads']} reads "
                  f"({result['pass_rate']:.1f}% pass)")


def example_process_with_map(sample_files: list[Path]):
    """Demonstrate process_with_map (memory-efficient)."""
    print("\n" + "=" * 60)
    print("Method 2: process_with_map (memory-efficient)")
    print("=" * 60)

    print("\nProcessing with map...")
    start_time = time.time()

    results = process_with_map(
        process_single_sample,
        sample_files,
        n_workers=3,
        chunksize=2  # Process 2 samples per worker at a time
    )

    elapsed = time.time() - start_time

    print(f"\nCompleted in {elapsed:.2f} seconds")
    print(f"Processed {len(results)} samples")


def example_parallel_map(sample_files: list[Path]):
    """Demonstrate parallel_map (simple interface)."""
    print("\n" + "=" * 60)
    print("Method 3: parallel_map (simple interface)")
    print("=" * 60)

    print("\nProcessing with parallel map...")
    start_time = time.time()

    results = parallel_map(
        process_single_sample,
        sample_files,
        n_workers=3
    )

    elapsed = time.time() - start_time

    print(f"\nCompleted in {elapsed:.2f} seconds")


def example_threads_vs_processes(sample_files: list[Path]):
    """Compare thread-based vs process-based parallelism."""
    print("\n" + "=" * 60)
    print("Threads vs Processes Comparison")
    print("=" * 60)

    # Process-based (better for CPU-bound tasks)
    print("\nUsing processes (CPU-bound tasks)...")
    start_time = time.time()
    results_processes = process_samples_parallel(
        sample_files,
        process_single_sample,
        n_workers=3,
        use_threads=False,
        show_progress=False
    )
    process_time = time.time() - start_time

    # Thread-based (better for I/O-bound tasks)
    print("Using threads (I/O-bound tasks)...")
    start_time = time.time()
    results_threads = process_samples_parallel(
        sample_files,
        process_single_sample,
        n_workers=3,
        use_threads=True,
        show_progress=False
    )
    thread_time = time.time() - start_time

    print(f"\nProcess-based: {process_time:.2f}s")
    print(f"Thread-based: {thread_time:.2f}s")


def process_with_potential_error(sample_file: Path) -> dict[str, int]:
    """Process function that might raise an error."""
    # Simulate error for specific sample
    if "sample_3" in sample_file.name:
        raise ValueError(f"Simulated error for {sample_file.name}")

    return {"sample": sample_file.name, "status": "success"}


def example_error_handling(sample_files: list[Path]):
    """Demonstrate error handling in parallel processing."""
    print("\n" + "=" * 60)
    print("Error Handling")
    print("=" * 60)

    print("\nProcessing with potential errors...")
    results = process_samples_parallel(
        sample_files,
        process_with_potential_error,
        n_workers=2,
        show_progress=False
    )

    print("\nResults:")
    for i, result in enumerate(results, 1):
        if result is None:
            print(f"  Sample {i}: ERROR (see above for details)")
        else:
            print(f"  Sample {i}: {result['status']}")


def main():
    """Run all parallel processing examples."""
    # Setup
    example_dir = Path("/tmp/parallel_example")
    print("Creating sample FASTQ files...")
    sample_files = create_sample_files(example_dir, n_samples=5)
    print(f"Created {len(sample_files)} sample files\n")

    # Run examples
    example_worker_detection()
    example_process_samples_parallel(sample_files)
    example_process_with_map(sample_files)
    example_parallel_map(sample_files)
    example_threads_vs_processes(sample_files)
    example_error_handling(sample_files)

    print("\n" + "=" * 60)
    print("Parallel Processing Best Practices:")
    print("=" * 60)
    print("1. Use processes (use_threads=False) for CPU-bound tasks")
    print("2. Use threads (use_threads=True) for I/O-bound tasks")
    print("3. Use process_samples_parallel for progress tracking")
    print("4. Use process_with_map for large datasets (more memory efficient)")
    print("5. Let get_optimal_workers() choose worker count automatically")
    print("6. Handle errors gracefully (results will be None for failed items)")
    print("=" * 60)


if __name__ == "__main__":
    main()
