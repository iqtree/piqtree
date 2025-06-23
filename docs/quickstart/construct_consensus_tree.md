# Construct a consensus tree from a collection of trees

A consensus tree can be constructed from a collection of `cogent3` tree objects with 
[`consensus_tree`](../api/tree/consensus_tree.md). The method defaults to the majority rule consensus 
tree - though different minimum levels of support can be used to achieve the strict consensus tree, 
and extended majority rule consensus tree.

## Usage

### Majority rule consensus tree

From a collection of trees, the majority rule consensus tree can be constructed.

```python
from cogent3 import make_tree
from piqtree import consensus_tree

tree1 = make_tree("(a,(b,(c,(d,e))));")
tree2 = make_tree("(a,((b,c),(d,e)));")
tree3 = make_tree("(a,((b,d),(c,e)));")

tree = consensus_tree([tree1, tree2, tree3]) # (a,(b,c,(d,e)));
```

### Strict consensus tree

The `min_support` parameter can be set to 1.0 construct a strict consensus tree (only the clades that appear in all of the trees).
The parameter represents the proportion of trees a clade must appear in to be in the resulting tree.

```python
from cogent3 import make_tree
from piqtree import consensus_tree

tree1 = make_tree("(a,(b,(c,(d,e))));")
tree2 = make_tree("(a,((b,c),(d,e)));")
tree3 = make_tree("(a,((b,d),(c,e)));")

tree = consensus_tree([tree1, tree2, tree3], min_support=1.0) # (a,(b,c,d,e));
```

### Extended majority rule consensus tree

When the `min_support` parameter is equal to 0, the extended majority rule consensus tree is computed. 

```python
from cogent3 import make_tree
from piqtree import consensus_tree

tree1 = make_tree("(a,(b,(c,(d,e))));")
tree2 = make_tree("(a,((b,c),(d,e)));")
tree3 = make_tree("(a,((b,c),(d,e)));")
tree4 = make_tree("(a,(c,(b,(d,e))));")
tree5 = make_tree("(c,(b,(a,(d,e))));")

tree = consensus_tree([tree1, tree2, tree3, tree4, tree5], min_support=0) # (a,((b,c),(d,e));
```

## See also

- For making a random phylogenetic tree, see ["Make a randomly generated phylogenetic tree"](make_random_tree.md).
