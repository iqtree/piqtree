# Fit branch lengths to a tree topology from an alignment

Branch lengths can be fitted to a tree from a cogent3 alignment object
using [`fit_tree`](../api/tree/fit_tree.md).

## Usage

### Basic Usage

Construct `cogent3` alignment and tree objects, then fit branch lengths to a new tree.

```python
from cogent3 import load_aligned_seqs, make_tree
from piqtree import fit_tree

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")
tree = make_tree("((Human, Chimpanzee), Rhesus, Mouse);")

fitted_tree = fit_tree(aln, tree, model="F81")
log_likelihood = fitted_tree.params["lnL"]
```

### Multithreading

To speed up computation, the number of threads to be used may be specified.
By default, the computation is done on a single thread. If 0 is specified,
then IQ-TREE attempts to determine the optimal number of threads.

> **Caution:** If 0 is specified with small datasets, the time to determine the
> optimal number of threads may exceed the time to find the maximum likelihood
> tree.

```python
from cogent3 import load_aligned_seqs, make_tree
from piqtree import fit_tree

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")
tree = make_tree("((Human, Chimpanzee), Rhesus, Mouse);")
model = "HKY+I+R3"

fitted_tree = fit_tree(aln, tree, model, num_threads=4)
```

### Fixed branch lengths

The likelihood for a tree with specified branch lengths may be calculated
by using the `bl_fixed` argument. When set to True, specified branch lengths
are not optimised and the likelihood for the tree is calculated. Any unspecified
branch lengths will still be optimised.

```python
from cogent3 import load_aligned_seqs, make_tree
from piqtree import fit_tree

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")
tree = make_tree("((Human:0.1, Chimpanzee:0.15):0.2, Rhesus:0.3, Mouse:0.5);")

tree_with_likelihood = fit_tree(aln, tree, model="GTR", bl_fixed=True)
log_likelihood = tree_with_likelihood.params["lnL"]
```

## See also

- For how to specify a `Model`, see ["Use different kinds of substitution models"](using_substitution_models.md).
- For constructing a maximum likelihood tree, see ["Construct a maximum likelihood phylogenetic tree"](construct_ml_tree.md).
