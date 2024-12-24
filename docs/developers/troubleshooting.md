# Troubleshooting

This is a guide for common problems that may occur when working on `piqtree` development. If you have an issue that isn't listed here,
please feel more than welcome to [raise an issue](https://github.com/iqtree/piqtree/issues) and we'll look into it.

## `piqtree` doesn't install

There could be several reasons `piqtree` doesn't install correctly.

### Check your operating system is supported

At this stage, we currently only support x86-64 linux, and x86-64 and ARM macOS.

### Check the IQ-TREE library is up-to-date

Try running the following to build the latest version of the IQ-TREE library.

```bash
git submodule update --init --recursive
./build_tools/build_iqtree.sh
```

If the above doesn't work, please re-try the guide on ["Setting up your environment for development"](./environment_setup.md).
