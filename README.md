# BioSeqFlow

> A comprehensive sequencing quality control and preprocessing platform

[![CI](https://github.com/omicpilot/bioseqflow/workflows/CI/badge.svg)](https://github.com/omicpilot/bioseqflow/actions)
[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

## Overview

BioSeqFlow is a user-friendly, modular platform designed to democratize access to sequencing quality control and preprocessing tools. Built for researchers with no coding background, it provides automated pipelines for various sequencing technologies including NGS, long-read sequencing, and more.

### Key Features

- **Multi-Platform Support**: NGS (Illumina), long-read (PacBio, Oxford Nanopore), and expandable to other sequencing technologies
- **Automated QC Pipeline**: FastQC, MultiQC, and custom quality metrics
- **Preprocessing Suite**: Adapter trimming, quality filtering, deduplication, and contamination screening
- **Parallel Processing**: Efficient batch processing of multiple samples
- **Interactive Visualizations**: Publication-ready plots with Plotly and Seaborn
- **User-Friendly CLI**: Simple command-line interface for non-programmers
- **Comprehensive Reports**: HTML reports with embedded interactive plots
- **Modular Design**: Easy to extend with custom modules

## Installation

### From PyPI (coming soon)

```bash
pip install bioseqflow
```

### From Source

```bash
git clone https://github.com/omicpilot/bioseqflow.git
cd bioseqflow
pip install -e ".[dev]"
```

### Development Setup

```bash
# Clone the repository
git clone https://github.com/omicpilot/bioseqflow.git
cd bioseqflow

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode with all dependencies
pip install -e ".[dev,dashboard,docs]"

# Install pre-commit hooks
pre-commit install
```

## Quick Start

### Basic QC Analysis

```bash
# Run quality control on a single FASTQ file
bioseqflow qc --input sample.fastq.gz --output results/

# Run with custom parameters
bioseqflow qc --input sample.fastq.gz --output results/ --threads 8 --min-quality 25
```

### Full Pipeline (QC + Preprocessing)

```bash
# Run complete pipeline
bioseqflow pipeline \
    --input sample.fastq.gz \
    --output results/ \
    --adapter AGATCGGAAGAG \
    --min-quality 20 \
    --min-length 50
```

### Batch Processing

```bash
# Process multiple samples from a sample sheet
bioseqflow batch --sample-sheet samples.csv --output results/
```

## Project Structure

```
bioseqflow/
├── bioseqflow/              # Main package
│   ├── core/               # Core functionality
│   │   ├── base.py        # Base classes for modules
│   │   └── config.py      # Configuration management
│   ├── qc/                # Quality control modules
│   │   ├── fastqc.py     # FastQC integration
│   │   ├── metrics.py    # Custom QC metrics
│   │   └── parsers.py    # Result parsers
│   ├── preprocessing/     # Preprocessing modules
│   │   ├── trimming.py   # Adapter trimming
│   │   ├── filtering.py  # Quality filtering
│   │   └── deduplication.py  # Duplicate removal
│   ├── utils/            # Utility functions
│   │   ├── io.py        # File I/O operations
│   │   ├── validators.py # Input validation
│   │   └── parallel.py   # Parallel processing
│   ├── visualization/    # Visualization modules
│   │   ├── plots.py     # Quality plots
│   │   └── reports.py   # HTML report generation
│   └── cli/             # Command-line interface
│       └── main.py      # CLI entry point
├── tests/               # Comprehensive test suite
├── docs/               # Documentation
├── examples/           # Example workflows and data
└── .github/workflows/  # CI/CD pipelines
```

## Usage Examples

### Python API

```python
from bioseqflow import Config
from bioseqflow.qc import FastQCRunner, QualityMetrics
from bioseqflow.preprocessing import AdapterTrimmer

# Configure the pipeline
config = Config(
    min_quality=20,
    min_length=50,
    threads=4
)

# Run quality control
qc_runner = FastQCRunner(config)
qc_results = qc_runner.run("sample.fastq.gz")

# Calculate custom metrics
metrics = QualityMetrics()
quality_score = metrics.calculate_composite_score(qc_results)

# Trim adapters
trimmer = AdapterTrimmer(config)
trimmer.trim("sample.fastq.gz", "trimmed.fastq.gz")
```

### Sample Sheet Format

```csv
sample_id,fastq_r1,fastq_r2,adapter,min_quality
sample1,data/s1_R1.fastq.gz,data/s1_R2.fastq.gz,AGATCGGAAGAG,20
sample2,data/s2_R1.fastq.gz,data/s2_R2.fastq.gz,AGATCGGAAGAG,20
```

## Pipeline Stages

### 1. Initial QC
- Run FastQC on raw reads
- Parse quality metrics programmatically
- Generate custom quality statistics
- Identify potential issues

### 2. Preprocessing
- Adapter trimming (Trimmomatic/fastp)
- Quality filtering (Phred score-based)
- Duplicate removal (exact and clustering-based)
- Length filtering
- Contamination screening

### 3. Post-QC
- Rerun FastQC on processed reads
- Compare before/after statistics
- Generate MultiQC report
- Create interactive visualizations

### 4. Reporting
- HTML reports with embedded plots
- PDF summaries
- Tab-delimited statistics
- Interactive dashboard (optional)

## Biological Context

Understanding sequencing quality is critical because:

- **Poor quality bases** lead to false variant calls and incorrect biological conclusions
- **Adapter contamination** creates alignment artifacts and phantom sequences
- **PCR duplicates** introduce coverage bias and inflate variant allele frequencies
- **Proper QC** guides downstream analysis decisions and ensures reproducible results

## Advanced Features

### Composite Quality Score

BioSeqFlow implements a custom composite quality score that combines multiple metrics:

```python
def calculate_composite_quality_score(fastqc_data, trimming_stats):
    """
    Calculate a weighted composite quality score.

    Combines:
    - Per-base sequence quality
    - Per-sequence quality scores
    - Sequence length distribution
    - GC content deviation
    - Adapter content
    - Duplication levels

    Returns a score from 0-100, where >90 is excellent, 70-90 is good,
    50-70 needs review, and <50 requires attention.
    """
    # Implementation in bioseqflow.qc.metrics
```

### Parallel Processing

```python
from bioseqflow.utils import process_samples_parallel

# Process multiple samples in parallel
results = process_samples_parallel(
    sample_list=["s1.fastq", "s2.fastq", "s3.fastq"],
    function=run_qc,
    n_workers=4
)
```

## Development

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=bioseqflow --cov-report=html

# Run specific test categories
pytest -m unit           # Unit tests only
pytest -m integration    # Integration tests only
pytest -m "not slow"     # Skip slow tests
```

### Code Quality

```bash
# Format code
ruff format bioseqflow tests

# Lint code
ruff check bioseqflow tests

# Type check
mypy bioseqflow

# Run all checks
pre-commit run --all-files
```

### Adding New Modules

1. Create your module in the appropriate directory
2. Inherit from `QCModule` or `PreprocessingModule` base classes
3. Implement required methods
4. Add comprehensive tests
5. Update documentation

Example:

```python
from bioseqflow.core.base import QCModule

class CustomQCModule(QCModule):
    """Custom quality control module."""

    def run(self, input_file: str) -> dict:
        """Run the QC analysis."""
        # Implementation
        pass
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes with tests
4. Run the test suite and linters
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## Roadmap

- [x] Core QC functionality
- [x] Basic preprocessing pipeline
- [ ] Long-read sequencing support (PacBio, ONT)
- [ ] Interactive dashboard (Streamlit/Dash)
- [ ] Cloud integration (AWS, GCP)
- [ ] RNA-seq specific QC modules
- [ ] Single-cell sequencing support
- [ ] Real-time processing for MinION
- [ ] Machine learning-based quality prediction

## Citation

If you use BioSeqFlow in your research, please cite:

```bibtex
@software{bioseqflow2024,
  title = {BioSeqFlow: A Comprehensive Sequencing Quality Control Platform},
  author = {OmicsPilot},
  year = {2024},
  url = {https://github.com/omicpilot/bioseqflow}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Built as part of a project-based bioinformatics learning journey
- Integrates industry-standard tools: FastQC, Trimmomatic, MultiQC
- Inspired by the need to democratize bioinformatics for wet-lab researchers

## Support

- **Documentation**: [Coming soon]
- **Issues**: [GitHub Issues](https://github.com/omicpilot/bioseqflow/issues)
- **Discussions**: [GitHub Discussions](https://github.com/omicpilot/bioseqflow/discussions)

---

Made by [OmicsPilot](https://github.com/omicpilot)
