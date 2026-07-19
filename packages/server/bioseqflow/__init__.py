from __future__ import annotations

"""
BioSeqFlow: A comprehensive sequencing quality control and preprocessing platform.

This package provides tools for quality control, preprocessing, and visualization
of sequencing data from various platforms (NGS, long-read, etc.).
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("bioseqflow")
except PackageNotFoundError:
    __version__ = "unknown"

__author__ = "OmicsPilot"

from bioseqflow.core.base import PreprocessingModule, QCModule
from bioseqflow.core.config import Config

__all__ = [
    "Config",
    "QCModule",
    "PreprocessingModule",
    "__version__",
]
