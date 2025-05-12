# Release checklist

These instructions are for `piqtree` release managers.
The release is largely automated by the CI, but there are several things to confirm and do in this process.

The release process is automated through the "Release" GitHub Action.
The documentation is fetched by readthedocs from the "Build Docs" GitHub Action.

- The `piqtree` version has been correctly bumped.
- The testing, linting and type checking all pass on all supported platforms.
- The code has been thoroughly tested (check what's been missed in the coverage report).
- The documentation builds and appears correct.
- The documentation has been updated on readthedocs (this must be triggered from readthedocs).
- The "Release" GitHub Action has correctly uploaded to Test PyPI.
- The "Release" GitHub action has correctly uploaded to PyPI.
