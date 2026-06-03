# Construct a consensus tree from a collection of trees

A consensus tree can be constructed from a collection of `cogent3` tree objects with 
[`consensus_tree`](../api/tree/consensus_tree.md). The method defaults to the majority-rule consensus 
tree - though different minimum levels of support can be used to achieve the strict consensus tree, 
and extended majority-rule consensus tree. By default, the support values of the input trees are 
carried over onto the consensus tree (see [Bootstrap values](#bootstrap-values) below).

## Usage

### Majority-rule consensus tree

From a collection of trees, the majority-rule consensus tree can be constructed.

```python
from cogent3 import make_tree
from piqtree import consensus_tree

tree1 = make_tree("(a,(b,(c,(d,e))));")
tree2 = make_tree("(a,((b,c),(d,e)));")
tree3 = make_tree("(a,((b,d),(c,e)));")

tree = consensus_tree([tree1, tree2, tree3]) # (a,(b,c,(d,e)));
```

### Strict consensus tree

The `min_support` parameter can be set to 1.0 to construct a strict consensus tree (only the clades that appear in all of the trees).
The parameter represents the proportion of trees a clade must appear in to be in the resulting tree.

```python
from cogent3 import make_tree
from piqtree import consensus_tree

tree1 = make_tree("(a,(b,(c,(d,e))));")
tree2 = make_tree("(a,((b,c),(d,e)));")
tree3 = make_tree("(a,((b,d),(c,e)));")

tree = consensus_tree([tree1, tree2, tree3], min_support=1.0) # (a,(b,c,d,e));
```

### Extended majority-rule consensus tree

When the `min_support` parameter is equal to 0, the extended majority-rule consensus tree is computed. 

```python
from cogent3 import make_tree
from piqtree import consensus_tree

tree1 = make_tree("(a,(b,(c,(d,e))));")
tree2 = make_tree("(a,((b,c),(d,e)));")
tree3 = make_tree("(a,((b,c),(d,e)));")
tree4 = make_tree("(a,(c,(b,(d,e))));")
tree5 = make_tree("(c,(b,(a,(d,e))));")

tree = consensus_tree([tree1, tree2, tree3, tree4, tree5], min_support=0) # (a,((b,c),(d,e)));
```

### Bootstrap values

When the input trees carry their own branch support - for example, the ultrafast bootstrap support
produced by [`build_tree`](construct_ml_tree.md) - `consensus_tree` carries these values over onto
the matching clades of the consensus tree, averaged across the input trees that contain each clade.
The way the values are combined is set by `support_aggregate` (`"mean"` by default, or `"max"` /
`"min"`), and each node's value can be accessed from `#!py3 node.support`.

```python
from cogent3 import make_tree
from piqtree import consensus_tree

# the numbers on the internal nodes are each tree's own bootstrap support
tree1 = make_tree("(a,((b,c)96,(d,e)90));")
tree2 = make_tree("(a,((b,c)88,(d,e)70));")

tree = consensus_tree([tree1, tree2])  # support_aggregate="mean" by default
# (b,c) -> mean(96, 88) -> 92.0
# (d,e) -> mean(90, 70) -> 80.0
```

## See also

- For making a random phylogenetic tree, see ["Make a randomly generated phylogenetic tree"](make_random_tree.md).
