# Construct a maximum likelihood phylogenetic tree

A maximum-likelihood phylogenetic tree can be constructed from
cogent3 alignment objects using [`build_tree`](../api/tree/build_tree.md).

## Usage

### Basic Usage

Construct a `cogent3` alignment object, then construct a maximum-likelihood tree.

```python
from cogent3 import load_aligned_seqs
from piqtree import build_tree

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

tree = build_tree(aln, "GTR")
log_likelihood = tree.params["lnL"]
```

### Ultrafast Bootstrap

To perform ultrafast bootstrapping, the number of replicates can be specified. The number of replicates must be at least 1000.
The support for each node in the tree object can be accessed from `#!py3 node.params["support"]`.

```python
from cogent3 import load_aligned_seqs
from piqtree import build_tree

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")
tree = build_tree(aln, "K81+FO", bootstrap_replicates=2000)
```

### Reproducible Results

For reproducible results, a random seed may be specified.
> **Caution:** 0 and None are equivalent to no random seed being specified.

```python
from cogent3 import load_aligned_seqs
from piqtree import build_tree

aln = load_aligned_seqs("my_protein.fasta", moltype="protein")

tree = build_tree(aln, "Dayhoff", rand_seed=3)
```

### Multithreading

To speed up computation, the number of threads to be used may be specified.
By default, the computation is done on a single thread. If 0 is specified,
then IQ-TREE attempts to determine the optimal number of threads.

> **Caution:** If 0 is specified with small datasets, the time to determine the
> optimal number of threads may exceed the time to find the maximum likelihood
> tree.

```python
from cogent3 import load_aligned_seqs
from piqtree import build_tree

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")
model = "HKY+I+G6"

tree = build_tree(aln, model, num_threads=4)
```

## See also

- For how to specify a `Model`, see ["Use different kinds of substitution models"](using_substitution_models.md).
- For selecting the best `Model`, see ["Find the model of best fit with ModelFinder"](using_model_finder.md).
- For fitting branch lengths to a tree topology see ["Fit branch lengths to a tree topology from an alignment"](fit_tree_topology.md).
