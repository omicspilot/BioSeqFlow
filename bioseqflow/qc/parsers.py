from __future__ import annotations

"""Parsers for QC tool outputs."""

from pathlib import Path
from typing import Any


class FastQCParser:
    """Parser for FastQC output files."""

    def parse(self, fastqc_data_path: Path | str) -> dict[str, Any]:
        """
        Parse FastQC data file.

        Args:
            fastqc_data_path: Path to fastqc_data.txt file

        Returns:
            Dictionary with parsed data
        """
        fastqc_data_path = Path(fastqc_data_path)

        if not fastqc_data_path.exists():
            raise FileNotFoundError(f"FastQC data file not found: {fastqc_data_path}")

        results: dict[str, Any] = {
            "basic_statistics": {},
            "modules": {},
            "summary": {},
        }

        current_module = None
        module_data: list[str] = []

        with open(fastqc_data_path) as f:
            for line in f:
                line = line.strip()

                if line.startswith(">>"):
                    # Save previous module data
                    if current_module and module_data:
                        results["modules"][current_module] = self._parse_module_data(
                            current_module, module_data
                        )
                        module_data = []

                    # Start new module
                    parts = line[2:].split("\t")
                    current_module = parts[0]

                    if len(parts) > 1:
                        results["summary"][current_module] = parts[1]

                elif line.startswith("#"):
                    # Skip comments
                    continue

                elif current_module:
                    # Collect module data
                    module_data.append(line)

        # Save last module
        if current_module and module_data:
            results["modules"][current_module] = self._parse_module_data(
                current_module, module_data
            )

        return results

    def _parse_module_data(
        self, module_name: str, data_lines: list[str]
    ) -> dict[str, Any]:
        """
        Parse individual module data.

        Args:
            module_name: Name of the module
            data_lines: Lines of data for this module

        Returns:
            Parsed module data
        """
        if module_name == "Basic Statistics":
            return self._parse_basic_statistics(data_lines)

        # For other modules, store as list of lines for now
        # Can be expanded with specific parsers for each module
        return {"raw_data": data_lines}

    def _parse_basic_statistics(self, data_lines: list[str]) -> dict[str, str | int]:
        """
        Parse Basic Statistics module.

        Args:
            data_lines: Lines from Basic Statistics section

        Returns:
            Dictionary of basic statistics
        """
        stats: dict[str, str | int] = {}

        for line in data_lines:
            if "\t" in line:
                key, value = line.split("\t", 1)
                # Try to convert to int if possible
                try:
                    stats[key] = int(value)
                except ValueError:
                    stats[key] = value

        return stats


class MultiQCParser:
    """Parser for MultiQC output files."""

    def parse(self, multiqc_data_dir: Path | str) -> dict[str, Any]:
        """
        Parse MultiQC data directory.

        Args:
            multiqc_data_dir: Path to multiqc_data directory

        Returns:
            Dictionary with parsed MultiQC data
        """
        multiqc_data_dir = Path(multiqc_data_dir)

        if not multiqc_data_dir.exists():
            raise FileNotFoundError(f"MultiQC data directory not found: {multiqc_data_dir}")

        results: dict[str, Any] = {
            "general_stats": {},
            "fastqc_data": {},
        }

        # Look for general stats file
        general_stats_file = multiqc_data_dir / "multiqc_general_stats.txt"
        if general_stats_file.exists():
            results["general_stats"] = self._parse_general_stats(general_stats_file)

        # Look for FastQC-specific data
        fastqc_file = multiqc_data_dir / "multiqc_fastqc.txt"
        if fastqc_file.exists():
            results["fastqc_data"] = self._parse_fastqc_data(fastqc_file)

        return results

    def _parse_general_stats(self, stats_file: Path) -> dict[str, Any]:
        """
        Parse MultiQC general statistics file.

        Args:
            stats_file: Path to general stats file

        Returns:
            Dictionary of general statistics
        """
        stats: dict[str, Any] = {}

        with open(stats_file) as f:
            header = f.readline().strip().split("\t")

            for line in f:
                values = line.strip().split("\t")
                sample_name = values[0]
                stats[sample_name] = dict(zip(header[1:], values[1:]))

        return stats

    def _parse_fastqc_data(self, fastqc_file: Path) -> dict[str, Any]:
        """
        Parse MultiQC FastQC data file.

        Args:
            fastqc_file: Path to FastQC data file

        Returns:
            Dictionary of FastQC data
        """
        data: dict[str, Any] = {}

        with open(fastqc_file) as f:
            header = f.readline().strip().split("\t")

            for line in f:
                values = line.strip().split("\t")
                sample_name = values[0]
                data[sample_name] = dict(zip(header[1:], values[1:]))

        return data
