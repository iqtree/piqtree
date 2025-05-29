# piqtree documentation

## Overview

`piqtree` is a Python package that exposes selected [IQ-TREE](http://www.iqtree.org) capabilities within Python, using the [cogent3](https://cogent3.org) library as the interface.

`piqtree` is implemented with the goals of:

- making the individual high-performance components of IQ-TREE available within Python, enabling the community to take advantage of these routines.
- facilitating exploratory analyses by leveraging cogent3's capabilities to provide a rich user experience in interactive Jupyter notebooks, e.g. trivial parallelisation across collections of alignments.
- code using piqtree apps should be easy to understand.

In addition to the functions provided, `piqtree` provides mini-applications in the form of [cogent3 apps](https://cogent3.org/doc/app/index.html). These can interplay with other such apps, e.g. the [cogent3-ete3](https://pypi.org/project/cogent3-ete3/) tree conversion plugin, the [diverse-seqs](https://pypi.org/project/diverse-seq/) sequence subsampling plugin.

> **Note**
> `piqtree` does not implement all of the capabilities of IQ-TREE!

## Installation

You get the vanilla version of `piqtree` by running the following command.

```bash
pip install piqtree
```

To get visualisation support with plotly, use the `[extra]` option.

```bash
pip install "piqtree[extra]"
```

## Supported platforms

At present we do not provide native binaries for Windows. Windows users can run `piqtree` using the Windows Subsystem for Linux (WSL) which can installed via the Windows Store.
