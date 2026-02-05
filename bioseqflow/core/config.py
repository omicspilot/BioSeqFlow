from __future__ import annotations

"""Configuration management for BioSeqFlow."""

from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, Field, field_validator


class Config(BaseModel):
    """Main configuration class for BioSeqFlow."""

    # Quality filtering parameters
    min_quality: int = Field(default=20, ge=0, le=93, description="Minimum Phred quality score")
    min_length: int = Field(default=50, ge=1, description="Minimum read length after trimming")
    max_length: int | None = Field(default=None, description="Maximum read length")

    # Adapter parameters
    adapter_sequence: str | None = Field(
        default=None, description="Adapter sequence for trimming"
    )
    adapter_r2: str | None = Field(
        default=None, description="Adapter sequence for R2 (paired-end)"
    )

    # Processing parameters
    threads: int = Field(default=1, ge=1, description="Number of threads to use")
    chunk_size: int = Field(default=10000, ge=1000, description="Chunk size for processing")

    # Output parameters
    output_dir: Path = Field(default=Path("./output"), description="Output directory")
    keep_intermediate: bool = Field(
        default=False, description="Keep intermediate files"
    )
    compress_output: bool = Field(default=True, description="Compress output files")

    # QC thresholds
    min_base_quality: int = Field(default=20, ge=0, le=93)
    gc_content_min: float = Field(default=0.2, ge=0.0, le=1.0)
    gc_content_max: float = Field(default=0.8, ge=0.0, le=1.0)
    duplication_threshold: float = Field(default=0.7, ge=0.0, le=1.0)

    class Config:
        """Pydantic configuration."""

        arbitrary_types_allowed = True

    @field_validator("output_dir")
    @classmethod
    def validate_output_dir(cls, v: Path) -> Path:
        """Validate and create output directory if needed."""
        if isinstance(v, str):
            v = Path(v)
        return v

    @classmethod
    def from_yaml(cls, yaml_path: Path | str) -> "Config":
        """
        Load configuration from YAML file.

        Args:
            yaml_path: Path to YAML configuration file

        Returns:
            Config object

        Raises:
            FileNotFoundError: If YAML file doesn't exist
            ValueError: If YAML is invalid
        """
        yaml_path = Path(yaml_path)
        if not yaml_path.exists():
            raise FileNotFoundError(f"Config file not found: {yaml_path}")

        with open(yaml_path) as f:
            config_dict = yaml.safe_load(f)

        return cls(**config_dict)

    def to_yaml(self, yaml_path: Path | str) -> None:
        """
        Save configuration to YAML file.

        Args:
            yaml_path: Path to save YAML file
        """
        yaml_path = Path(yaml_path)
        yaml_path.parent.mkdir(parents=True, exist_ok=True)

        with open(yaml_path, "w") as f:
            yaml.safe_dump(
                self.model_dump(mode="python"),
                f,
                default_flow_style=False,
                sort_keys=False,
            )

    def to_dict(self) -> dict[str, Any]:
        """
        Convert configuration to dictionary.

        Returns:
            Configuration as dictionary
        """
        return self.model_dump(mode="python")

    def update(self, **kwargs: Any) -> None:
        """
        Update configuration parameters.

        Args:
            **kwargs: Parameters to update
        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
