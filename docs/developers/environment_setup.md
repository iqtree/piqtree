# Setting up your environment for development

We currently only support x86-64 linux, and both x86-64 and ARM mac operating systems.
For windows development we recommend using Windows Subsytem for Linux. Alternatively,
one can work through the provided dev container.

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

If you are working on linux, check the top of the `.devconatiner/Dockerfile` for the list of
dependencies. Once installed, the IQ-TREE library can be build with `./build_tools/build_iqtree.sh`.

## Installing `piqtree` for standard development

After completing the above steps, in your preferred Python virtual environment run the following
in the repository's directory:

```bash
pip install -e ".[dev]"
```

This will install `piqtree` in editable mode, alongside testing, linting and other dependencies.

## Installing `piqtree` for documentation development

To contribute to the documentation, we recommend running the following after the above
steps in a separate virtual environment:

```bash
pip install -e ".[doc]"
```

This will install `piqtree` in editable mode, along with dependencies for building the documentation.

## Running the tests

To verify that installation has worked, using your chosen standard development environment run
the following in the base directory of this repository.

```bash
pytest
```
