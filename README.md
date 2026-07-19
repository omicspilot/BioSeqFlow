<div align="center">
  <h1>BioSeqFlow</h1>
  <h2>NGS Quality Control & Preprocessing Pipeline</h2>
  <p>
    <b>A comprehensive, type-safe platform for automated sequencing quality control and preprocessing. Trim adapters, filter reads, detect duplicates, and analyze paired-end data — with fuzzy matching and parallel processing built in.</b>
  </p>
  <p>
    <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT"></a>
    <a href="https://pypi.org/project/bioseqflow/"><img src="https://img.shields.io/pypi/v/bioseqflow.svg" alt="PyPI version"></a>
    <a href="https://codecov.io/gh/omicspilot/BioSeqFlow"><img src="https://codecov.io/gh/omicspilot/BioSeqFlow/graph/badge.svg" alt="Codecov"></a>
    <br/>
    <a href="https://github.com/omicspilot/BioSeqFlow/actions/workflows/ci.yml"><img src="https://github.com/omicspilot/BioSeqFlow/actions/workflows/ci.yml/badge.svg" alt="CI workflow"></a>
    <a href="https://github.com/omicspilot/BioSeqFlow/commits/main"><img src="https://img.shields.io/github/last-commit/omicspilot/BioSeqFlow/main" alt="GitHub last commit"></a>
  </p>
</div>

---

## Overview

BioSeqFlow is a modular, production-ready platform for automated sequencing quality control and preprocessing. Built with type safety and comprehensive testing (74% coverage, 176 tests), it provides researchers with reliable tools for NGS data analysis.

### Key Features

- **Paired-End Support** - Complete R1/R2 validation, insert size analysis, and orientation checking
- **Fuzzy Adapter Matching** - Smith-Waterman alignment detects 30-50% more adapters than exact matching
- **QC Integration** - Programmatic FastQC/MultiQC parsing and analysis
- **Parallel Processing** - Multiple strategies with automatic worker optimization
- **Type-Safe** - Full mypy validation with comprehensive type hints
- **Well-Tested** - 176 tests covering edge cases and error conditions
- **Modular Design** - Easy to extend with custom preprocessing modules

## Installation

```bash
# From source (recommended for development)
git clone https://github.com/omicpilot/bioseqflow.git
cd bioseqflow/packages/server
pip install -e ".[dev]"

# From PyPI (coming soon)
pip install bioseqflow
```

## Quick Start

```python
from bioseqflow.utils.paired_end import validate_paired_files
from bioseqflow.preprocessing.trimming import AdapterTrimmer

# Validate paired-end files
result = validate_paired_files("sample_R1.fastq.gz", "sample_R2.fastq.gz")
print(f"OK: {result.total_pairs} valid pairs" if result.is_valid else f"ERROR: {result.errors}")

# Fuzzy adapter trimming (30-50% better detection)
trimmer = AdapterTrimmer()
stats = trimmer.trim(
    "sample.fastq.gz",
    "trimmed.fastq.gz",
    adapter="AGATCGGAAGAGC",
    fuzzy_match=True,        # Enable Smith-Waterman
    max_error_rate=0.15      # Tolerate 15% mismatches
)
```

## Documentation & Examples

**Full Documentation**: [omicspilot.com/projects/bioseqflow](https://omicspilot.com/projects/bioseqflow)

**API Reference**: [docs/usage.md](docs/usage.md)

## Project Structure

This repo is a monorepo: the QC/preprocessing package below is being grown into
**BioSeqFlow Desktop**, a GUI companion app that wraps the full NGS toolchain
(FastQC, fastp, STAR, BWA-MEM2, samtools, GATK, featureCounts, pyDESeq2, SnpEff)
behind a FastAPI server and an Electron + React UI. See
[context/CONTEXT.md](context/CONTEXT.md) and [context/ROADMAP.md](context/ROADMAP.md)
for the full architecture and implementation plan.

```
bioseqflow/
├── packages/
│   ├── server/                        # pip-installable Python package
│   │   ├── bioseqflow/                # QC/preprocessing toolkit (this README's features)
│   │   │   ├── core/                 # Base classes and configuration
│   │   │   ├── qc/                   # QC modules (FastQC, MultiQC, parsers)
│   │   │   ├── preprocessing/        # Trimming, filtering, deduplication
│   │   │   ├── utils/                # Paired-end, parallel, alignment utilities
│   │   │   ├── visualization/        # Plots and reports
│   │   │   └── cli/                  # Command-line interface
│   │   ├── bioseqflow_server/         # FastAPI app scaffold (subprocess orchestration, WIP)
│   │   ├── tests/                     # 176 comprehensive tests (74% coverage)
│   │   ├── pyproject.toml
│   │   └── pytest.ini
│   └── desktop/                       # Electron + React GUI scaffold (WIP)
│       ├── src/
│       └── electron/
└── docs/                              # API documentation and guides
```

## Core Capabilities

### Paired-End Sequencing
- R1/R2 concordance validation
- Multiple FASTQ header format support (Illumina, SRA)
- Insert size distribution calculation
- Read orientation consistency checking

### Fuzzy Adapter Matching
- Smith-Waterman local alignment
- Tolerates mismatches and partial matches
- Configurable error rates
- 30-50% improvement over exact matching

### QC & Analysis
- FastQC/MultiQC integration
- Programmatic results parsing
- Custom quality metrics
- Composite quality scoring

### Parallel Processing
- Automatic worker detection
- Process and thread-based strategies
- Progress tracking and error handling
- Memory-efficient batch processing

## Roadmap

### In Progress
- Interactive dashboard (Streamlit/Dash)
- MultiQC report generation
- Additional visualization modules

### Planned
- Long-read sequencing support (PacBio, ONT)
- Cloud integration (AWS, GCP)
- RNA-seq specific QC modules
- Single-cell sequencing support

## Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass and coverage is maintained
5. Submit a pull request

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

## Citation

If you use BioSeqFlow in your research, please cite:

```bibtex
@software{bioseqflow2025,
  title = {BioSeqFlow: A Comprehensive Sequencing Quality Control Platform},
  author = {OmicsPilot},
  year = {2024},
  url = {https://github.com/omicpilot/bioseqflow}
}
```

## Support

- **Discussions**: [GitHub Discussions](https://github.com/omicpilot/bioseqflow/discussions)
- **Issues**: [GitHub Issues](https://github.com/omicpilot/bioseqflow/issues)
