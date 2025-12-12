# Make a randomly generated phylogenetic tree

A randomly generated phylogenetic tree can be made using [`random_tree`](../api/tree/random_tree.md#piqtree.random_tree).
Multiple tree generation modes are supported with the [`TreeGenMode`](../api/tree/random_tree.md#piqtree.TreeGenMode) including
balanced, birth-death, caterpillar, star, uniform, and Yule-Harding tree.

## Usage

### Basic Usage

Specify the number of taxa, and under what method the tree is to be generated.
See the documentation for [`TreeGenMode`](../api/tree/random_tree.md#piqtree.TreeGenMode) for all available generation options.

> **Note:** if star trees are generated the tree appears bifurcating, but branch lengths are set to zero where appropriate.

```python
from piqtree import TreeGenMode, random_tree

num_taxa = 100

tree = random_tree(num_taxa, TreeGenMode.YULE_HARDING)
```

### Reproducible Results

For reproducible results, a random seed may be specified.
> **Caution:** 0 is a specific random seed. None is equivalent to no random seed being specified.

```python
from piqtree import TreeGenMode, random_tree

num_taxa = 100

tree = random_tree(num_taxa, TreeGenMode.CATERPILLAR, rand_seed=1)
```

## See also

- For constructing a maximum likelihood tree, see ["Construct a maximum likelihood phylogenetic tree"](construct_ml_tree.md).
- For constructing a consensus tree, see ["Construct a consensus tree from a collection of trees"](construct_consensus_tree.md).