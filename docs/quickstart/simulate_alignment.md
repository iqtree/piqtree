# Simulate alignments with AliSim

An alignment can be simulated over a cogent3 tree object using [`simulate_alignment`](../api/alignment/simulate_alignment.md).

## Usage

### Basic Usage

Construct a `cogent3` tree object, then simulate an alignment over it.

```python
from cogent3 import make_tree
from piqtree import simulate_alignment

tree = make_tree("(A:0.3544,(B:0.1905,C:0.1328):0.0998,D:0.0898);")
aln = simulate_alignment(tree, "JC")
```

### Reproducible Results

For reproducible results, a random seed may be specified.
> **Caution:** 0 is a specific random seed. None is equivalent to no random seed being specified.

```python
from cogent3 import make_tree
from piqtree import simulate_alignment

tree = make_tree("(A:0.3544,(B:0.1905,C:0.1328):0.0998,D:0.0898);")
aln = simulate_alignment(tree, "Dayhoff", rand_seed=42)
```

### Multithreading

To speed up computation, the number of threads to be used may be specified.
By default, the computation is done on a single thread. If 0 is specified,
then IQ-TREE automatically determines the optimal number of threads.

```python
from cogent3 import make_tree
from piqtree import simulate_alignment

tree = make_tree("(A:0.3544,(B:0.1905,C:0.1328):0.0998,D:0.0898);")
aln = simulate_alignment(tree, "JC", length=10000, num_threads=0)
```

## See also

- For how to specify a `Model`, see ["Use different kinds of substitution models"](using_substitution_models.md).
- For constructing a maximum likelihood tree, see ["Construct a maximum likelihood phylogenetic tree"](construct_ml_tree.md).