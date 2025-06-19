# Fit branch lengths to a tree topology from an alignment

Branch lengths can be fitted to a tree from a cogent3 alignment object
using [`fit_tree`](../api/tree/fit_tree.md).

## Usage

### Basic Usage

Construct `cogent3` alignment and tree objects, then fit branch lengths to a new tree.

```python
from cogent3 import load_aligned_seqs, make_tree
from piqtree import Model, fit_tree

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")
tree = make_tree("((Human, Chimpanzee), Rhesus, Mouse);")

fitted_tree = fit_tree(aln, tree, model=Model("F81"))
log_likelihood = fitted_tree.params["lnL"]
```

### Reproducible Results

For reproducible results, a random seed may be specified.
> **Caution:** 0 and None are equivalent to no random seed being specified.

```python
from cogent3 import load_aligned_seqs, make_tree
from piqtree import Model, fit_tree
from piqtree.model import StandardDnaModel

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")
tree = make_tree("((Human, Chimpanzee), Rhesus, Mouse);")

fitted_tree = fit_tree(aln, tree, model=Model(StandardDnaModel.SYM), rand_seed=42)
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
from piqtree import Model, fit_tree
from piqtree.model import StandardDnaModel, FreeRateModel

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")
tree = make_tree("((Human, Chimpanzee), Rhesus, Mouse);")

model = Model(StandardDnaModel.HKY, rate_model=FreeRateModel(3), invariant_sites=True)
fitted_tree = fit_tree(aln, tree, model, num_threads=4)
```

## See also

- For how to specify a `Model`, see ["Use different kinds of substitution models"](using_substitution_models.md).
- For constructing a maximum likelihood tree, see ["Construct a maximum likelihood phylogenetic tree"](construct_ml_tree.md).
