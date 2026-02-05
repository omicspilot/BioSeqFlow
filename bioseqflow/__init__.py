from __future__ import annotations

"""
BioSeqFlow: A comprehensive sequencing quality control and preprocessing platform.

This package provides tools for quality control, preprocessing, and visualization
of sequencing data from various platforms (NGS, long-read, etc.).
"""

__version__ = "0.1.0"
__author__ = "OmicsPilot"

from bioseqflow.core.config import Config
from bioseqflow.core.base import QCModule, PreprocessingModule

__all__ = [
    "Config",
    "QCModule",
    "PreprocessingModule",
    "__version__",
]
