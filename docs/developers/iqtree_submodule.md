# Working with the IQ-TREE submodule

We use a submodule to keep track of the version of IQ-TREE being used.

Every so often, there is an update to the submodule.
This most frequently happens when new features are being added from IQ-TREE, or
upon a new release.

When the submodule updates, to ensure the latest version of the IQ-TREE library is being used
the submodule must be updated locally and the library rebuilt. This can be done by running the
following two commands in the base directory of the `piqtree` repository:

```bash
git submodule update --init --recursive
./build_tools/build_iqtree.sh
```
