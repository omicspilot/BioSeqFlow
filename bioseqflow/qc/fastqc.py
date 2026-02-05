from __future__ import annotations

"""FastQC integration module."""

import subprocess
from pathlib import Path
from typing import Any

from bioseqflow.core.base import QCModule
from bioseqflow.qc.parsers import FastQCParser
from bioseqflow.utils.validators import validate_file_path, validate_threads


class FastQCRunner(QCModule):
    """Run FastQC analysis on sequencing files."""

    def __init__(self, config: Any | None = None) -> None:
        """
        Initialize FastQC runner.

        Args:
            config: Configuration object
        """
        super().__init__(config)
        self.parser = FastQCParser()

    def validate_inputs(self, input_file: Path | str, output_dir: Path | str) -> None:
        """
        Validate input parameters.

        Args:
            input_file: Input FASTQ file
            output_dir: Output directory

        Raises:
            FileNotFoundError: If input file doesn't exist
            ValueError: If inputs are invalid
        """
        validate_file_path(
            input_file,
            must_exist=True,
            extensions=[".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        )

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    def run(
        self,
        input_file: Path | str,
        output_dir: Path | str,
        threads: int = 1,
        quiet: bool = False,
    ) -> dict[str, Any]:
        """
        Run FastQC on input file.

        Args:
            input_file: Input FASTQ file
            output_dir: Output directory for FastQC results
            threads: Number of threads to use
            quiet: Suppress FastQC output

        Returns:
            Dictionary containing FastQC results and statistics

        Raises:
            RuntimeError: If FastQC execution fails
        """
        self.validate_inputs(input_file, output_dir)
        validate_threads(threads)

        input_file = Path(input_file)
        output_dir = Path(output_dir)

        cmd = [
            "fastqc",
            str(input_file),
            "--outdir",
            str(output_dir),
            "--threads",
            str(threads),
        ]

        if quiet:
            cmd.append("--quiet")

        try:
            result = subprocess.run(
                cmd, check=True, capture_output=True, text=True
            )

            # Parse the output
            self.results = self.parse_results(output_dir)
            self.results["fastqc_stdout"] = result.stdout

            return self.results

        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"FastQC failed with error: {e.stderr}"
            ) from e
        except FileNotFoundError as e:
            raise RuntimeError(
                "FastQC not found. Please install FastQC and ensure it's in your PATH."
            ) from e

    def parse_results(self, output_path: Path) -> dict[str, Any]:
        """
        Parse FastQC output files.

        Args:
            output_path: Directory containing FastQC output

        Returns:
            Parsed results dictionary
        """
        # Find the fastqc_data.txt file
        data_files = list(output_path.glob("**/fastqc_data.txt"))

        if not data_files:
            return {"error": "FastQC data file not found"}

        return self.parser.parse(data_files[0])


def run_fastqc(
    input_file: Path | str,
    output_dir: Path | str,
    threads: int = 1,
    quiet: bool = True,
) -> dict[str, Any]:
    """
    Convenience function to run FastQC.

    Args:
        input_file: Input FASTQ file
        output_dir: Output directory
        threads: Number of threads
        quiet: Suppress output

    Returns:
        FastQC results dictionary
    """
    runner = FastQCRunner()
    return runner.run(input_file, output_dir, threads=threads, quiet=quiet)
