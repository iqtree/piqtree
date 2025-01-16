# piqtree

[![PyPI Version](https://img.shields.io/pypi/v/piqtree)](https://pypi.org/project/piqtree/)
[![Python Version](https://img.shields.io/pypi/pyversions/piqtree)](https://pypi.org/project/piqtree/)
[![License](https://img.shields.io/github/license/iqtree/piqtree)](https://github.com/iqtree/piqtree/blob/main/LICENSE)

[![CI](https://github.com/iqtree/piqtree/workflows/CI/badge.svg)](https://github.com/iqtree/piqtree/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/iqtree/piqtree/badge.svg?branch=main)](https://coveralls.io/github/iqtree/piqtree?branch=main)
[![Documentation Status](https://readthedocs.org/projects/piqtree/badge/?version=latest)](https://piqtree.readthedocs.io/en/latest/?badge=latest)

`piqtree` is a library which allows you use IQ-TREE directly from Python! The interface with python is through [cogent3](https://cogent3.org) objects.
For usage, please refer to the [documentation](https://piqtree.readthedocs.io/) or the examples below.

If you encounter any problems or have any feature requests feel free to raise an [issue](https://github.com/iqtree/piqtree/issues)!

## Contributing

If you would like to help out by contributing to the piqtree project, please check out our [contributor guide!](https://piqtree.readthedocs.io/en/latest/developers/)

## Examples

### Phylogenetic Reconstruction

```python
from piqtree import build_tree
from piqtree.model import Model
from cogent3 import load_aligned_seqs # Included with piqtree!

# Load Sequences
aln = load_aligned_seqs("tests/data/example.fasta", moltype="dna")
aln = aln.take_seqs(["Human", "Chimpanzee", "Rhesus", "Mouse"])

# Reconstruct a phylogenetic tree with IQ-TREE!
tree = build_tree(aln, Model("JC"), rand_seed=1) # Optionally specify a random seed.

print("Tree topology:", tree) # A cogent3 tree object
print("Log-likelihood:", tree.params["lnL"])
# In a Jupyter notebook, try tree.get_figure() to see a dendrogram
```

> **Note**
> See the [cogent3 docs](https://cogent3.org) for examples on what you can do with cogent3 trees.

### Fit Branch Lengths to Tree Topology

```python
from piqtree import fit_tree
from piqtree.model import Model
from cogent3 import load_aligned_seqs, make_tree # Included with piqtree!

# Load Sequences
aln = load_aligned_seqs("tests/data/example.fasta", moltype="dna")
aln = aln.take_seqs(["Human", "Chimpanzee", "Rhesus", "Mouse"])

# Construct tree topology
tree = make_tree("(Human, Chimpanzee, (Rhesus, Mouse));")

# Fit branch lengths with IQ-TREE!
tree = fit_tree(aln, tree, Model("JC"), rand_seed=1) # Optionally specify a random seed.

print("Tree with branch lengths:", tree) # A cogent3 tree object
print("Log-likelihood:", tree.params["lnL"])
```

### Create a Collection of Random Trees

```python
from piqtree import TreeGenMode, random_trees

num_trees = 3 
num_taxa = 5

# Also supports YULE_HARDING, CATERPILLAR, BALANCED, BIRTH_DEATH and STAR_TREE
tree_gen_mode = TreeGenMode.UNIFORM 

# Randomly generate trees
trees = random_trees(num_trees, num_taxa, tree_gen_mode, rand_seed=1) # Optionally specify a random seed.

print(trees) # A tuple of 3 trees with 5 taxa each.
```

### Pairwise Robinson-Foulds Distance between Trees

```python
from piqtree import robinson_foulds
from cogent3 import make_tree # Included with piqtree!

# Construct trees
tree1 = make_tree("(a,b,(c,(d,e)));")
tree2 = make_tree("(e,b,(c,(d,a)));")
tree3 = make_tree("(a,b,(d,(c,e)));")

# Calculate pairwise distances
pairwise_distances = robinson_foulds([tree1, tree2, tree3]) # Supports any number of trees (for a sequence of trees use *seq_of_trees)

print(pairwise_distances) # A numpy array containing pairwaise Robinson-Foulds distances between trees
```
