.PHONY: dev check test clean

# First-time setup: install dev dependencies and git hooks
dev:
	uv sync --extra dev
	pre-commit install
	pre-commit install --hook-type commit-msg

# Run the full quality pass (mirrors the pre-commit hooks + CI lint job)
check:
	uv run mypy bioseqflow_server
	uv run ruff check --fix bioseqflow_server tests
	uv run ruff format bioseqflow_server tests
	uv run black bioseqflow_server tests

# Run the test suite
test:
	uv run pytest

# Remove the virtual environment (re-run `make dev` to rebuild)
clean:
	rm -rf .venv
