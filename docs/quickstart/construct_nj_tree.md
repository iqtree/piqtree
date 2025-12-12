# Construct a rapid neighbour-joining tree from a distance matrix

A rapid neighbour-joining tree can be constructed from a pairwise `DistanceMatrix` object with [`nj_tree`](../api/tree/nj_tree.md).

## Usage

### Basic Usage

Make a pairwise distance matrix, then construct the rapid neighbour-joining tree.

```python
from cogent3 import load_aligned_seqs
from piqtree import jc_distances, nj_tree

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

distance_matrix = jc_distances(aln)

tree = nj_tree(distance_matrix)
```

### Other Distance Matrices

`cogent3` supports the calculation of **paralinear**, **JC69**, **TN93**, **hamming** and **pdist** distance matrices from alignment objects.

```python
from cogent3 import load_aligned_seqs
from piqtree import nj_tree

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

distance_matrix = aln.distance_matrix(calc="TN93") # calc is not case sensitive

tree = nj_tree(distance_matrix)
```

## See also

- For constructing the JC distance matrix with `piqtree`, see ["Construct a rapid neighbour-joining tree from a distance matrix"](calculate_jc_distances.md).
- For constructing a maximum likelihood tree, see ["Construct a maximum likelihood phylogenetic tree"](construct_ml_tree.md).
