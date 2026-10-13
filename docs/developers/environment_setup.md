# Setting up your environment for development

## Initial setup

Fork and clone the [`piqtree` repository](https://github.com/iqtree/piqtree).

## Working with the `iqtree3` submodule

We use a git submodule to keep track of the version of IQ-TREE that is being used.
To initialise the submodule, or to get the latest version after an update, the
following can be run:

```bash
git submodule update --init --recursive
```

## Building the IQ-TREE library

There are several build scripts used by the CI to install dependencies and build the
IQ-TREE library.

If you are working through the dev container, all dependencies should already be installed.
The IQ-TREE library can be built with `./build_tools/build_iqtree.sh`.

If you are working on a Mac, running `./build_tools/before_all_mac.sh` will install the
required dependencies through homebrew, then build the IQ-TREE library. If you need to
build the library again, running `./build_tools/build_iqtree.sh` will skip the dependency
installation step.

If you are working on linux, check the top of the `.devcontainer/DockerFile` for the list of
dependencies. Once installed, the IQ-TREE library can be built with `./build_tools/build_iqtree.sh`.

## Installing `piqtree` for standard development

We use [`uv`](https://docs.astral.sh/uv/) to manage the development environment. See the
[`uv` installation guide](https://docs.astral.sh/uv/getting-started/installation/) if you don't
already have it.

After completing the above steps, run the following in the repository's directory:

```bash
uv sync
```

This creates a virtual environment in `.venv` and installs `piqtree` in editable mode, alongside
testing, linting and other dependencies (the `dev` dependency group). The exact versions are pinned
by `uv.lock`. `uv` will download a suitable version of Python if needed.

Commands can be run inside the environment with `uv run`, e.g. `uv run pytest`, or by activating
it with `source .venv/bin/activate` (`.venv\Scripts\activate` on Windows).

`piqtree`'s C++ extension is compiled when it is installed, so changes to it (or to the IQ-TREE
library) require a rebuild. `uv sync` and `uv run` will rebuild it when any of the files listed under
`cache-keys` in `pyproject.toml` have changed. Running commands directly from an activated environment
does not check for changes, so run `uv sync` first. To force a rebuild, run:

```bash
uv sync --reinstall-package piqtree
```

## Installing `piqtree` for documentation development

To contribute to the documentation, also install the `doc` dependency group:

```bash
uv sync --group doc
```

This will install `piqtree` in editable mode, along with dependencies for building the documentation.

## Managing dependencies

Dependencies are declared in `pyproject.toml` and locked in `uv.lock`, which is committed to the repository.
Use `uv add` and `uv remove` to change dependencies (e.g. `uv add --group test <package>` for a testing dependency),
which update both files. If you edit `pyproject.toml` by hand, run `uv lock` and commit the updated `uv.lock`.
The CI will fail if `uv.lock` is out of date.

## Running the tests

To verify that installation has worked, run the following in the base directory of this repository.

```bash
uv run pytest
```
