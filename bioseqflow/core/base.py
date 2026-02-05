from __future__ import annotations

"""Base classes for BioSeqFlow modules."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any


class BaseModule(ABC):
    """Abstract base class for all BioSeqFlow modules."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize the module.

        Args:
            config: Configuration object or dictionary
        """
        self.config = config

    @abstractmethod
    def validate_inputs(self, *args: Any, **kwargs: Any) -> None:
        """
        Validate input parameters.

        Raises:
            ValueError: If inputs are invalid
        """
        pass

    @abstractmethod
    def run(self, *args: Any, **kwargs: Any) -> Any:
        """
        Execute the module's main functionality.

        Returns:
            Module-specific output
        """
        pass


class QCModule(BaseModule):
    """Base class for quality control modules."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize QC module.

        Args:
            config: Configuration object
        """
        super().__init__(config)
        self.results: dict[str, Any] = {}

    @abstractmethod
    def parse_results(self, output_path: Path) -> dict[str, Any]:
        """
        Parse the module's output files.

        Args:
            output_path: Path to output files

        Returns:
            Parsed results as dictionary
        """
        pass

    def get_summary(self) -> dict[str, Any]:
        """
        Get summary statistics from results.

        Returns:
            Summary statistics dictionary
        """
        return self.results.get("summary", {})


class PreprocessingModule(BaseModule):
    """Base class for preprocessing modules."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize preprocessing module.

        Args:
            config: Configuration object
        """
        super().__init__(config)
        self.stats: dict[str, int | float] = {}

    @abstractmethod
    def process(self, input_file: Path, output_file: Path) -> dict[str, int | float]:
        """
        Process the input file and write to output.

        Args:
            input_file: Input file path
            output_file: Output file path

        Returns:
            Processing statistics
        """
        pass

    def get_stats(self) -> dict[str, int | float]:
        """
        Get processing statistics.

        Returns:
            Statistics dictionary with counts and percentages
        """
        return self.stats
