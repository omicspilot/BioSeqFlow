# Contributing to BioSeqFlow

Thank you for your interest in contributing to BioSeqFlow! This document provides guidelines for contributing to the project.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/bioseqflow.git
   cd bioseqflow
   ```

3. **Set up development environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -e ".[dev]"
   pre-commit install
   ```

## Development Workflow

### 1. Create a Branch

Create a new branch for your feature or bugfix:

```bash
git checkout -b feature/your-feature-name
```

Branch naming conventions:
- `feature/` - New features
- `fix/` - Bug fixes
- `docs/` - Documentation improvements
- `test/` - Test additions or improvements
- `refactor/` - Code refactoring

### 2. Make Your Changes

Follow these guidelines:

#### Code Style

- **Follow PEP 8** guidelines
- **Use type hints** for all function signatures
- **Write docstrings** for all public functions and classes (Google style)
- **Keep functions small and focused** (max ~50 lines)
- **Use descriptive variable names**

Example:

```python
def calculate_quality_score(
    sequences: list[str],
    min_quality: int = 20
) -> dict[str, float]:
    """
    Calculate quality scores for sequences.

    Args:
        sequences: List of DNA sequences
        min_quality: Minimum quality threshold

    Returns:
        Dictionary with quality statistics

    Raises:
        ValueError: If sequences list is empty
    """
    if not sequences:
        raise ValueError("Sequences list cannot be empty")

    # Implementation here
    return {"mean_quality": 30.5}
```

#### Testing

- **Write tests** for all new code
- **Aim for >80% coverage**
- **Use pytest fixtures** for reusable test data
- **Test edge cases** and error conditions

Example:

```python
import pytest
from bioseqflow.utils.validators import validate_sequence

def test_validate_sequence():
    """Test sequence validation."""
    # Valid case
    validate_sequence("ATCGATCG")

    # Invalid case
    with pytest.raises(ValueError):
        validate_sequence("INVALID")
```

#### Documentation

- **Update docstrings** when changing function behavior
- **Add examples** to documentation
- **Update README.md** for significant features
- **Add type hints** to all function signatures

### 3. Run Tests and Checks

Before committing, ensure all checks pass:

```bash
# Format code
ruff format bioseqflow tests

# Lint code
ruff check bioseqflow tests

# Type check
mypy bioseqflow

# Run tests
pytest

# Run tests with coverage
pytest --cov=bioseqflow --cov-report=html

# Run all pre-commit hooks
pre-commit run --all-files
```

### 4. Commit Your Changes

Write clear, descriptive commit messages:

```bash
git add .
git commit -m "feat: add support for long-read sequencing

- Add PacBio adapter trimming
- Update quality score calculations
- Add tests for new functionality"
```

Commit message format:
- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation changes
- `test:` - Test additions/changes
- `refactor:` - Code refactoring
- `style:` - Code style changes
- `chore:` - Maintenance tasks

### 5. Push and Create Pull Request

```bash
git push origin feature/your-feature-name
```

Then create a Pull Request on GitHub with:
- **Clear title** describing the change
- **Description** of what changed and why
- **Link to related issues** (if applicable)
- **Screenshots** (for UI changes)
- **Breaking changes** noted clearly

## Code Review Process

1. Maintainers will review your PR
2. Address any feedback or requested changes
3. Once approved, your PR will be merged

## Testing Guidelines

### Unit Tests

Test individual functions and classes:

```python
def test_adapter_trimmer():
    """Test adapter trimming functionality."""
    trimmer = AdapterTrimmer()
    stats = trimmer.trim("input.fastq", "output.fastq", "AGATCGGAAGAG")

    assert stats["total_reads"] > 0
    assert stats["percent_trimmed"] >= 0
```

### Integration Tests

Test multiple components working together:

```python
@pytest.mark.integration
def test_complete_pipeline():
    """Test full QC pipeline."""
    # Run multiple steps
    # Verify end-to-end functionality
```

### Test Fixtures

Use fixtures for reusable test data:

```python
@pytest.fixture
def sample_config():
    """Sample configuration for testing."""
    return Config(min_quality=20, min_length=50)
```

## Adding New Features

### 1. Discuss First

For major features, open an issue first to discuss:
- Design approach
- API changes
- Breaking changes
- Implementation plan

### 2. Follow Module Structure

New modules should follow the existing pattern:

```python
from bioseqflow.core.base import QCModule

class NewQCModule(QCModule):
    """Description of the module."""

    def __init__(self, config=None):
        """Initialize module."""
        super().__init__(config)

    def validate_inputs(self, *args, **kwargs):
        """Validate input parameters."""
        pass

    def run(self, *args, **kwargs):
        """Execute the module."""
        pass

    def parse_results(self, output_path):
        """Parse output files."""
        pass
```

### 3. Add Tests

Include comprehensive tests:
- Unit tests for new functions
- Integration tests for workflows
- Edge case tests
- Error handling tests

### 4. Update Documentation

- Add docstrings
- Update usage guide
- Add examples
- Update README if needed

## Reporting Issues

### Bug Reports

Include:
- **Clear title and description**
- **Steps to reproduce**
- **Expected vs actual behavior**
- **Environment details** (OS, Python version, etc.)
- **Minimal reproducible example**
- **Error messages and stack traces**

### Feature Requests

Include:
- **Clear use case description**
- **Proposed solution**
- **Alternative approaches considered**
- **Potential impact** on existing functionality

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Accept constructive criticism
- Focus on what's best for the community
- Show empathy towards others

## Questions?

- Open an issue for questions
- Check existing issues and PRs
- Read the documentation

## License

By contributing, you agree that your contributions will be licensed under the same MIT License that covers the project.

---

Thank you for contributing to BioSeqFlow!
