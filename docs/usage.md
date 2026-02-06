# BioSeqFlow Usage Guide

## Installation

```bash
pip install -e ".[dev]"
```

## Quick Start

### 1. Basic QC Analysis

Run quality control on a FASTQ file:

```bash
bioseqflow qc --input sample.fastq.gz --output qc_results/ --threads 4
```

### 2. Adapter Trimming

Trim adapter sequences:

```bash
bioseqflow trim \
    --input sample.fastq.gz \
    --output trimmed.fastq.gz \
    --adapter AGATCGGAAGAG
```

### 3. Quality Filtering

Filter reads by quality and length:

```bash
bioseqflow filter \
    --input sample.fastq.gz \
    --output filtered.fastq.gz \
    --min-quality 20 \
    --min-length 50
```

### 4. Complete Pipeline

Run the full QC and preprocessing pipeline:

```bash
bioseqflow pipeline \
    --input sample.fastq.gz \
    --output results/ \
    --adapter AGATCGGAAGAG \
    --min-quality 20 \
    --min-length 50 \
    --threads 4
```

### 5. Batch Processing

Process multiple samples from a sample sheet:

```bash
bioseqflow batch \
    --sample-sheet samples.csv \
    --output batch_results/ \
    --threads 8
```

## Python API

### Basic Usage

```python
from bioseqflow import Config
from bioseqflow.qc import FastQCRunner
from bioseqflow.preprocessing import AdapterTrimmer, QualityFilter

# Create configuration
config = Config(
    min_quality=20,
    min_length=50,
    threads=4,
    adapter_sequence="AGATCGGAAGAG"
)

# Run QC
qc_runner = FastQCRunner(config)
results = qc_runner.run("sample.fastq.gz", "qc_output/")

# Trim adapters
trimmer = AdapterTrimmer(config)
stats = trimmer.trim("sample.fastq.gz", "trimmed.fastq.gz", "AGATCGGAAGAG")

# Filter reads
filter_module = QualityFilter(config)
filter_stats = filter_module.filter("trimmed.fastq.gz", "filtered.fastq.gz")
```

### Working with FASTQ Files

```python
from bioseqflow.utils.io import read_fastq, write_fastq

# Read FASTQ file
for record in read_fastq("sample.fastq.gz"):
    print(f"Name: {record.name}")
    print(f"Sequence: {record.sequence}")
    print(f"Length: {record.length}")
    print(f"Mean Quality: {record.mean_quality()}")

# Write FASTQ file
records = list(read_fastq("input.fastq"))
write_fastq(records, "output.fastq.gz", compress=True)
```

### Custom Quality Metrics

```python
from bioseqflow.qc.metrics import QualityMetrics

metrics = QualityMetrics()

# Calculate GC content
gc_content = metrics.calculate_gc_content("ATCGATCGATCG")

# Calculate composite quality score
score = metrics.calculate_composite_score(fastqc_results)
print(f"Quality Score: {score}/100")
```

### Parallel Processing

```python
from bioseqflow.utils.parallel import process_samples_parallel

def process_sample(sample_file):
    # Your processing function
    return run_qc(sample_file)

samples = ["s1.fastq", "s2.fastq", "s3.fastq"]
results = process_samples_parallel(
    samples,
    process_sample,
    n_workers=4,
    show_progress=True
)
```

### Paired-End Sequencing

Validate and process paired-end FASTQ files:

```python
from bioseqflow.utils.paired_end import (
    validate_paired_files,
    read_paired_fastq,
    calculate_insert_size_distribution,
    check_read_orientation
)

# Validate R1/R2 file pairing
result = validate_paired_files("sample_R1.fastq.gz", "sample_R2.fastq.gz")
if result.is_valid:
    print(f"✓ Files are properly paired ({result.total_pairs} pairs)")
else:
    print(f"✗ Validation errors: {result.errors}")

# Read paired FASTQ files simultaneously
for pair in read_paired_fastq("sample_R1.fastq.gz", "sample_R2.fastq.gz"):
    print(f"R1: {pair.r1.sequence}")
    print(f"R2: {pair.r2.sequence}")

# Calculate insert size distribution
stats = calculate_insert_size_distribution(
    "sample_R1.fastq.gz",
    "sample_R2.fastq.gz",
    max_reads=10000
)
print(f"Mean insert size: {stats['mean_insert_size']:.1f} bp")
print(f"Insert size range: {stats['min_insert_size']}-{stats['max_insert_size']} bp")

# Check read orientation consistency
orientation = check_read_orientation("sample_R1.fastq.gz", "sample_R2.fastq.gz")
if orientation['length_consistent']:
    print("✓ Read lengths are consistent")
else:
    print(f"⚠ Mean length difference: {orientation['mean_length_diff']:.1f} bp")
```

### Fuzzy Adapter Matching

Use Smith-Waterman alignment for more sensitive adapter detection:

```python
from bioseqflow.preprocessing.trimming import AdapterTrimmer
from bioseqflow.utils.alignment import smith_waterman, find_adapter_fuzzy

# Trim adapters with fuzzy matching (tolerates mismatches and partial matches)
trimmer = AdapterTrimmer()
stats = trimmer.trim(
    "sample.fastq.gz",
    "trimmed.fastq.gz",
    adapter="AGATCGGAAGAGC",
    fuzzy_match=True,       # Enable fuzzy matching
    max_error_rate=0.15     # Allow up to 15% mismatches
)
print(f"Trimmed {stats['trimmed_reads']} reads ({stats['trim_rate']:.1f}%)")

# Find adapter in a specific sequence
sequence = "ATCGATCGATCGATCGAGATCGGAAGAGCGTCGT"
adapter = "AGATCGGAAGAGC"

# Exact search
result = find_adapter_fuzzy(sequence, adapter, max_error_rate=0.0)

# Fuzzy search (allows mismatches)
result_fuzzy = find_adapter_fuzzy(sequence, adapter, max_error_rate=0.15)
if result_fuzzy:
    start, end = result_fuzzy
    print(f"Adapter found at position {start}-{end}: {sequence[start:end]}")

# Direct Smith-Waterman alignment
alignment = smith_waterman(adapter, sequence)
print(f"Alignment score: {alignment.score}")
print(f"Aligned adapter: {alignment.aligned_query}")
print(f"Aligned read:    {alignment.aligned_target}")
```

### Parsing QC Results

Parse FastQC and MultiQC output files:

```python
from bioseqflow.qc.parsers import FastQCParser, MultiQCParser

# Parse FastQC output
fastqc_parser = FastQCParser()
results = fastqc_parser.parse("qc_results/sample_fastqc/fastqc_data.txt")

# Access basic statistics
basic_stats = results['modules']['Basic Statistics']
print(f"Total sequences: {basic_stats['Total Sequences']}")
print(f"GC content: {basic_stats['%GC']}%")
print(f"Sequence length: {basic_stats['Sequence length']}")

# Check module status
for module, status in results['summary'].items():
    print(f"{module}: {status}")

# Parse MultiQC output
multiqc_parser = MultiQCParser()
multiqc_results = multiqc_parser.parse("multiqc_data/")

# Access per-sample statistics
for sample, stats in multiqc_results['general_stats'].items():
    print(f"{sample}: {stats}")
```

## Configuration

### YAML Configuration

Create a configuration file (`config.yaml`):

```yaml
min_quality: 20
min_length: 50
adapter_sequence: "AGATCGGAAGAG"
threads: 4
output_dir: "./results"
compress_output: true
```

Load it in Python:

```python
from bioseqflow.core.config import Config

config = Config.from_yaml("config.yaml")
```

### Sample Sheet Format

Sample sheet CSV format:

```csv
sample_id,fastq_r1,fastq_r2,adapter,min_quality
sample1,s1_R1.fastq.gz,s1_R2.fastq.gz,AGATCGGAAGAG,20
sample2,s2_R1.fastq.gz,s2_R2.fastq.gz,AGATCGGAAGAG,25
```

## Advanced Features

### Advanced Parallel Processing

BioSeqFlow provides multiple parallel processing strategies:

```python
from bioseqflow.utils.parallel import (
    process_samples_parallel,
    process_with_map,
    parallel_map,
    get_optimal_workers
)

# Method 1: process_samples_parallel (with progress bar and error handling)
def analyze_sample(fastq_file):
    # Your analysis function
    return run_analysis(fastq_file)

samples = ["s1.fastq", "s2.fastq", "s3.fastq"]

# Automatic worker detection
workers = get_optimal_workers()  # Uses CPU count - 1
print(f"Using {workers} workers")

# Process with custom workers and progress
results = process_samples_parallel(
    samples,
    analyze_sample,
    n_workers=4,
    use_threads=False,      # Use processes (CPU-bound) vs threads (I/O-bound)
    show_progress=True,
    description="Analyzing samples"
)

# Method 2: process_with_map (memory-efficient for large datasets)
results = process_with_map(
    analyze_sample,
    samples,
    n_workers=4,
    chunksize=2  # Process 2 items per worker at a time
)

# Method 3: parallel_map (simple parallel map)
results = parallel_map(analyze_sample, samples, n_workers=4)
```

### Deduplication

Remove duplicate reads:

```bash
bioseqflow deduplicate \
    --input sample.fastq.gz \
    --output dedup.fastq.gz \
    --method exact
```

In Python:

```python
from bioseqflow.preprocessing import DuplicateRemover

dedup = DuplicateRemover()
stats = dedup.remove_duplicates(
    "sample.fastq.gz",
    "dedup.fastq.gz",
    method="exact",
    keep_best=True
)
```

### Custom Workflows

```python
from pathlib import Path
from bioseqflow import Config
from bioseqflow.qc import FastQCRunner
from bioseqflow.preprocessing import (
    AdapterTrimmer,
    QualityTrimmer,
    QualityFilter,
    DuplicateRemover
)

# Setup
config = Config(min_quality=25, min_length=75, threads=8)
output_dir = Path("custom_workflow")
output_dir.mkdir(exist_ok=True)

# Step 1: QC
qc = FastQCRunner(config)
initial_qc = qc.run("input.fastq.gz", output_dir / "qc_initial")

# Step 2: Trim adapters
trimmer = AdapterTrimmer(config)
trimmer.trim("input.fastq.gz", output_dir / "trimmed.fastq.gz", "AGATCGGAAGAG")

# Step 3: Quality trim
qtrimmer = QualityTrimmer(config)
qtrimmer.trim(output_dir / "trimmed.fastq.gz", output_dir / "qtrimmed.fastq.gz")

# Step 4: Filter
filter_mod = QualityFilter(config)
filter_mod.filter(output_dir / "qtrimmed.fastq.gz", output_dir / "filtered.fastq.gz")

# Step 5: Deduplicate
dedup = DuplicateRemover(config)
dedup.remove_duplicates(output_dir / "filtered.fastq.gz", output_dir / "final.fastq.gz")

# Step 6: Final QC
final_qc = qc.run(output_dir / "final.fastq.gz", output_dir / "qc_final")
```

## Best Practices

1. **Always run QC before and after processing** to assess improvement
2. **Use appropriate quality thresholds** for your sequencing platform
3. **Keep intermediate files** during pipeline development
4. **Use parallel processing** for batch jobs
5. **Document your parameters** in configuration files
6. **Validate inputs** before processing large datasets

## Troubleshooting

### FastQC not found

Install FastQC and ensure it's in your PATH:

```bash
# macOS
brew install fastqc

# Linux
sudo apt-get install fastqc

# Or use conda
conda install -c bioconda fastqc
```

### Memory issues with large files

Use chunked processing or reduce the number of parallel workers:

```python
config = Config(chunk_size=5000, threads=2)
```

### Quality scores seem wrong

Check your quality score encoding (Phred+33 vs Phred+64):

```python
from bioseqflow.utils.validators import validate_quality_score

validate_quality_score(30, encoding="phred33")
```
