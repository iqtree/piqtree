# Testing, linting, and type checking

We use `pytest` for testing, `ruff` for linting and `mypy` for type checking. Pull requests require all three of these things to pass.

## Testing

To run the tests, run the following in the base directory of the `piqtree` repository:

```bash
pytest
```

## Linting

To run the linting, run the following in the base directory of the `piqtree` repository:

```bash
ruff check
```

## Type checking

To run type checking, run the following in the base directory of the `piqtree` repository:

```bash
mypy src tests
```
