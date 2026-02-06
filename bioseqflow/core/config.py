from __future__ import annotations

"""Configuration management for BioSeqFlow."""

from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator

# Organism-specific GC content profiles
GC_PROFILES = {
    "human": (0.35, 0.55),         # Homo sapiens: ~40-45% GC
    "mouse": (0.35, 0.50),         # Mus musculus: ~42% GC
    "bacteria": (0.25, 0.75),      # Broad range for bacteria
    "ecoli": (0.45, 0.55),         # E. coli: ~50% GC
    "yeast": (0.35, 0.45),         # S. cerevisiae: ~38% GC
    "arabidopsis": (0.32, 0.42),   # A. thaliana: ~36% GC
    "drosophila": (0.38, 0.48),    # D. melanogaster: ~42% GC
    "celegans": (0.32, 0.42),      # C. elegans: ~36% GC
    "malaria": (0.15, 0.25),       # P. falciparum: 19-20% GC (AT-rich)
    "tb": (0.60, 0.70),            # M. tuberculosis: 65% GC
    "streptomyces": (0.65, 0.75),  # Streptomyces: ~70% GC
    "fungi": (0.35, 0.60),         # Broad range for fungi
    "generic": (0.15, 0.90),       # Very permissive for unknown organisms
}


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
    gc_content_min: float = Field(
        default=0.15, ge=0.0, le=1.0, description="Minimum GC content (fraction)"
    )
    gc_content_max: float = Field(
        default=0.90, ge=0.0, le=1.0, description="Maximum GC content (fraction)"
    )
    duplication_threshold: float = Field(default=0.7, ge=0.0, le=1.0)
    organism: str | None = Field(
        default=None,
        description="Organism profile for GC thresholds (overrides gc_content_min/max)"
    )

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @field_validator("output_dir")
    @classmethod
    def validate_output_dir(cls, v: Path) -> Path:
        """Validate and create output directory if needed."""
        if isinstance(v, str):
            v = Path(v)
        return v

    def __init__(self, **data: Any) -> None:
        """Initialize config and apply organism-specific settings."""
        super().__init__(**data)

        # Apply organism-specific GC thresholds if organism is specified
        if self.organism:
            if self.organism in GC_PROFILES:
                gc_min, gc_max = GC_PROFILES[self.organism]
                self.gc_content_min = gc_min
                self.gc_content_max = gc_max
            else:
                raise ValueError(
                    f"Unknown organism: {self.organism}. "
                    f"Available profiles: {list(GC_PROFILES.keys())}"
                )

    @classmethod
    def from_yaml(cls, yaml_path: Path | str) -> Config:
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

        # Convert to dict and ensure Path objects are strings
        config_dict = self.model_dump(mode="python")
        # Convert Path objects to strings for YAML serialization
        if "output_dir" in config_dict and isinstance(config_dict["output_dir"], Path):
            config_dict["output_dir"] = str(config_dict["output_dir"])

        with open(yaml_path, "w") as f:
            yaml.safe_dump(
                config_dict,
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
