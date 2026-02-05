"""Tests for core modules."""

import tempfile
from pathlib import Path

import pytest
import yaml

from bioseqflow.core.base import PreprocessingModule, QCModule
from bioseqflow.core.config import Config


class TestConfig:
    """Tests for Config class."""

    def test_default_config(self) -> None:
        """Test default configuration."""
        config = Config()

        assert config.min_quality == 20
        assert config.min_length == 50
        assert config.threads == 1
        assert config.keep_intermediate is False

    def test_custom_config(self) -> None:
        """Test custom configuration."""
        config = Config(
            min_quality=30,
            min_length=100,
            threads=4,
            adapter_sequence="AGATCGGAAGAG",
        )

        assert config.min_quality == 30
        assert config.min_length == 100
        assert config.threads == 4
        assert config.adapter_sequence == "AGATCGGAAGAG"

    def test_config_validation(self) -> None:
        """Test configuration validation."""
        # Quality score out of range
        with pytest.raises(Exception):
            Config(min_quality=100)

        # Invalid min_length
        with pytest.raises(Exception):
            Config(min_length=0)

    def test_config_to_dict(self) -> None:
        """Test config to dictionary conversion."""
        config = Config(min_quality=25, threads=2)
        config_dict = config.to_dict()

        assert isinstance(config_dict, dict)
        assert config_dict["min_quality"] == 25
        assert config_dict["threads"] == 2

    def test_config_from_yaml(self, temp_dir: Path) -> None:
        """Test loading config from YAML."""
        yaml_file = temp_dir / "config.yaml"

        config_data = {
            "min_quality": 25,
            "min_length": 75,
            "threads": 8,
            "adapter_sequence": "AGATCGGAAGAG",
        }

        with open(yaml_file, "w") as f:
            yaml.safe_dump(config_data, f)

        config = Config.from_yaml(yaml_file)

        assert config.min_quality == 25
        assert config.min_length == 75
        assert config.threads == 8
        assert config.adapter_sequence == "AGATCGGAAGAG"

    def test_config_to_yaml(self, temp_dir: Path) -> None:
        """Test saving config to YAML."""
        config = Config(min_quality=30, threads=4)
        yaml_file = temp_dir / "config_output.yaml"

        config.to_yaml(yaml_file)

        assert yaml_file.exists()

        with open(yaml_file) as f:
            loaded_data = yaml.safe_load(f)

        assert loaded_data["min_quality"] == 30
        assert loaded_data["threads"] == 4

    def test_config_update(self) -> None:
        """Test updating configuration."""
        config = Config(min_quality=20)

        config.update(min_quality=30, threads=8)

        assert config.min_quality == 30
        assert config.threads == 8


class DummyQCModule(QCModule):
    """Dummy QC module for testing."""

    def validate_inputs(self, *args, **kwargs) -> None:  # type: ignore[no-untyped-def]
        """Validate inputs."""
        pass

    def run(self, *args, **kwargs) -> dict:  # type: ignore[no-untyped-def]
        """Run module."""
        return {"test": "results"}

    def parse_results(self, output_path: Path) -> dict:  # type: ignore[type-arg]
        """Parse results."""
        return {"parsed": True}


class DummyPreprocessingModule(PreprocessingModule):
    """Dummy preprocessing module for testing."""

    def validate_inputs(self, *args, **kwargs) -> None:  # type: ignore[no-untyped-def]
        """Validate inputs."""
        pass

    def process(self, input_file: Path, output_file: Path) -> dict[str, int | float]:
        """Process file."""
        return {"processed": 1}


class TestBaseClasses:
    """Tests for base classes."""

    def test_qc_module_initialization(self) -> None:
        """Test QC module initialization."""
        module = DummyQCModule()

        assert module.results == {}
        assert module.config is None

    def test_qc_module_with_config(self) -> None:
        """Test QC module with config."""
        config = Config(min_quality=25)
        module = DummyQCModule(config)

        assert module.config == config

    def test_qc_module_run(self) -> None:
        """Test QC module run method."""
        module = DummyQCModule()
        results = module.run()

        assert results == {"test": "results"}

    def test_preprocessing_module_initialization(self) -> None:
        """Test preprocessing module initialization."""
        module = DummyPreprocessingModule()

        assert module.stats == {}
        assert module.config is None

    def test_preprocessing_module_get_stats(self) -> None:
        """Test getting preprocessing stats."""
        module = DummyPreprocessingModule()
        module.stats = {"total_reads": 1000, "passed_reads": 950}

        stats = module.get_stats()

        assert stats["total_reads"] == 1000
        assert stats["passed_reads"] == 950
