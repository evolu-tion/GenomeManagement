# Changelog

All notable changes to this project will be documented in this file.

## [2.0.0] - 2026-01-18

### Added

- **PlantPAN3 Support**: Added `plantpan3` command for integrating with PlantPAN 3.0.
- **Testing**: Added `pytest` infrastructure with unit tests (`tests/test_seq_manage.py`) and CLI integration tests (`tests/test_cli.py`).
- **CI/CD**: Added GitHub Actions for automated testing (`python-package.yml`) and PyPI publishing (`publish.yml`).
- **Documentation**: Created comprehensive GitHub Wiki documentation in `docs/wiki`.
- **Packaging**: Standardized project structure (`src` layout) and configuration via `pyproject.toml`.

### Changed

- **Refactoring**: Completely refactored `seq_manage.py` to adhere to PEP 8 standards (e.g., `Fasta_manager` -> `FastaManager`).
- **API**: Updated all internal scripts (`promoter_retrieve.py`, `proteins_retrieve.py`, etc.) to use the new standardized API.
- **Entry Points**: Registered console scripts via `pyproject.toml` instead of standalone script execution.
- **Installation**: Updated installation methods to support `pip install .` and PyPI.
- **Versioning**: Bumped version to `2.0.0` across all files.

### Fixed

- **CI Compatibility**: Fixed CLI tests to run correctly in CI environments by invoking python modules directly.
- **Dependencies**: Explicitly defined dependencies and optional test dependencies.
