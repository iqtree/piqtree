# Testing, linting, and type checking

We use `pytest` for testing, `ruff` for linting and `mypy` for type checking. Pull requests require all three of these things to pass.

## Testing

To run the tests, run the following in the base directory of the `piqtree` repository:

```bash
uv run pytest
```

## Linting

To run the linting, run the following in the base directory of the `piqtree` repository:

```bash
uv run ruff check
```

## Type checking

To run type checking, run the following in the base directory of the `piqtree` repository:

```bash
uv run mypy src tests
```

## Running checks across Python versions

The CI runs the above checks through [`nox`](https://nox.thea.codes/), which uses `uv` to create a fresh
environment from `uv.lock` for each supported Python version. To list the available sessions, run:

```bash
uv run nox -l
```

A specific session can then be run with, e.g., `uv run nox -s test-3.14`.
