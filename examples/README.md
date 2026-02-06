# BioSeqFlow Examples

This directory contains comprehensive examples demonstrating BioSeqFlow's key features.

## Available Examples

### 1. Fuzzy Adapter Trimming (`fuzzy_adapter_example.py`)

Demonstrates fuzzy adapter matching using Smith-Waterman alignment:
- Basic Smith-Waterman alignment
- Fuzzy adapter finding with error tolerance
- Comparison of exact vs fuzzy trimming performance

**Key Learning**: Fuzzy matching catches 30-50% more adapters than exact matching!

```bash
python examples/fuzzy_adapter_example.py
```

### 2. Paired-End Sequencing (`paired_end_example.py`)

Complete workflow for paired-end data validation and processing:
- Read ID extraction from various FASTQ header formats (Illumina old/new, SRA)
- R1/R2 file validation and concordance checking
- Simultaneous paired reading
- Insert size distribution calculation
- Read orientation consistency checks

```bash
python examples/paired_end_example.py
```

**Features Demonstrated**:
- `validate_paired_files()` - Ensure R1/R2 concordance
- `read_paired_fastq()` - Read pairs simultaneously
- `calculate_insert_size_distribution()` - Estimate insert sizes
- `check_read_orientation()` - Verify read length consistency
- `extract_read_id()` - Parse various header formats

### 3. Parallel Processing (`parallel_processing_example.py`)

Advanced parallel processing strategies for batch analysis:
- Automatic worker detection based on CPU count
- Process-based vs thread-based parallelism
- Multiple processing strategies (with progress, map-based, simple)
- Error handling in parallel workflows
- Performance comparisons

```bash
python examples/parallel_processing_example.py
```

**Methods Demonstrated**:
- `process_samples_parallel()` - With progress bar and error handling
- `process_with_map()` - Memory-efficient for large datasets
- `parallel_map()` - Simple parallel map
- `get_optimal_workers()` - Automatic worker detection

### 4. Parsing QC Results (`parsing_qc_results_example.py`)

Parse and analyze FastQC and MultiQC output files:
- Parsing FastQC data files programmatically
- Extracting quality metrics and module statuses
- Parsing MultiQC aggregated results across samples
- Sample comparison (R1 vs R2)
- Quality-based sample filtering

```bash
python examples/parsing_qc_results_example.py
```

**Features Demonstrated**:
- `FastQCParser` - Parse individual FastQC outputs
- `MultiQCParser` - Parse aggregated MultiQC data
- Extracting basic statistics
- Module status checking (pass/warn/fail)
- Cross-sample quality comparisons
- Automated quality filtering

## Configuration Files

### `config_example.yaml`

Example YAML configuration file showing all available parameters:
- Quality thresholds
- Adapter sequences
- Thread settings
- Output options

```yaml
min_quality: 20
min_length: 50
adapter_sequence: "AGATCGGAAGAG"
threads: 4
output_dir: "./results"
compress_output: true
```

### `sample_sheet_example.csv`

Example sample sheet for batch processing:
```csv
sample_id,fastq_r1,fastq_r2,adapter,min_quality
sample1,s1_R1.fastq.gz,s1_R2.fastq.gz,AGATCGGAAGAG,20
sample2,s2_R1.fastq.gz,s2_R2.fastq.gz,AGATCGGAAGAG,25
```

## Sample Data

The `sample_data/` directory contains small test FASTQ files for experimentation:
- Single-end samples
- Paired-end samples (R1/R2)
- Various quality profiles

## Notebooks

The `notebooks/` directory contains Jupyter notebooks with interactive tutorials:
- Quality control workflows
- Preprocessing pipelines
- Data visualization
- Batch processing

## Running Examples

### Prerequisites

```bash
# Install BioSeqFlow with development dependencies
pip install -e ".[dev]"
```

### Run All Examples

```bash
# Run individual examples
python examples/fuzzy_adapter_example.py
python examples/paired_end_example.py
python examples/parallel_processing_example.py
python examples/parsing_qc_results_example.py
```

### Using in Your Own Code

All examples are self-contained and can be adapted for your own workflows:

```python
# Copy and modify example code
from bioseqflow.utils.paired_end import validate_paired_files

# Use with your own data
result = validate_paired_files("my_R1.fastq.gz", "my_R2.fastq.gz")
if result.is_valid:
    print(f"✓ Ready to process {result.total_pairs} pairs")
```

## Best Practices Demonstrated

1. **Type Safety**: All examples use proper type hints
2. **Error Handling**: Graceful handling of edge cases
3. **Progress Tracking**: User-friendly progress bars for long operations
4. **Documentation**: Clear docstrings and comments
5. **Reproducibility**: Configuration files for reproducible workflows
6. **Performance**: Efficient parallel processing strategies

## Learning Path

Recommended order for learning:

1. **Start with**: `fuzzy_adapter_example.py` - Learn about alignment algorithms
2. **Then**: `paired_end_example.py` - Understand paired-end data
3. **Next**: `parsing_qc_results_example.py` - Learn to analyze QC outputs
4. **Finally**: `parallel_processing_example.py` - Scale up with parallel processing

## Getting Help

- Check the main documentation: `docs/usage.md`
- Review test files in `tests/` for more examples
- Open an issue on GitHub for specific questions

## Contributing Examples

Have a useful workflow? Contribute an example:

1. Create a self-contained Python script
2. Add clear docstrings and comments
3. Include example output or expected results
4. Update this README with your example
5. Submit a pull request!

---

**Need more help?** See the full documentation in the `docs/` directory or run:

```bash
python -c "from bioseqflow import __version__; print(__version__)"
```
