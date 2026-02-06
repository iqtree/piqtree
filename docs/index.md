## Overview

`piqtree` (pronounced 'pie-cue-tree') is a Python package that exposes selected [IQ-TREE](http://www.iqtree.org) capabilities within Python, using the [cogent3](https://cogent3.org) library as the interface.

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

To get Jupyter and visualisation support (with plotly) use the `[extra]` option.

```bash
pip install "piqtree[extra]"
```

## Sharing your cool solutions 😎

If you've done something you would like to share which either used `piqtree` or helped you to it (e.g. some data sampling steps with `cogent3`), please share your code on [our community contribution site](https://github.com/cogent3/c3codeshare). Contributing is as easy as clicking a button and pasting your code! Who knows, it may lead to collaborations 🙂.

If you want to see how others are using `piqtree`, take a look at the site. (At the time of writing, this is very new, so there's not many contributions.)

## Getting help

If you have a question about using `piqtree`, you can post a question [in the piqtree forums](https://github.com/iqtree/piqtree/discussions). Likewise, if you have any questions about using `cogent3`, post a question [in the cogent3 forums](https://github.com/cogent3/cogent3/discussions).

## Reporting problems

To report problems or potential bugs 🐛, please raise an issue on the [piqtree GitHub issues page](https://github.com/iqtree/piqtree/issues).

## Citing piqtree

Please cite our [preprint](https://www.biorxiv.org/content/10.1101/2025.07.13.664626v3).
