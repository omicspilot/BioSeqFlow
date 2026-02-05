"""Pytest configuration and shared fixtures."""

from pathlib import Path
from typing import Generator
import tempfile

import pytest


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_fastq_content() -> str:
    """Return sample FASTQ content for testing."""
    return """@SEQ_ID_1
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
@SEQ_ID_2
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
"""


@pytest.fixture
def sample_fastq_file(temp_dir: Path, sample_fastq_content: str) -> Path:
    """Create a sample FASTQ file for testing."""
    fastq_file = temp_dir / "sample.fastq"
    fastq_file.write_text(sample_fastq_content)
    return fastq_file


@pytest.fixture
def sample_config() -> dict[str, str | int]:
    """Return sample configuration for testing."""
    return {
        "min_quality": 20,
        "min_length": 50,
        "adapter_sequence": "AGATCGGAAGAG",
        "threads": 4,
    }
