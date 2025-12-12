# Calculate pairwise Jukes-Cantor distances

A pairwise Jukes-Cantor distance matrix can be constructed from a cogent3 alignment object using [`jc_distances`](../api/genetic_distance/jc_distances.md).
The resulting distance matrix can be indexed by integer index, or by name.

## Usage

### Basic Usage

Construct a `cogent3` alignment object, then calculate the pairwise JC distance matrix.

```python
from cogent3 import load_aligned_seqs
from piqtree import jc_distances

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

distance_matrix = jc_distances(aln)

distance_1 = distance_matrix[0, 1]

distance_2 = distance_matrix["Human", "Chimpanzee"]
```

### Multithreading

The number of threads to be used may be specified. By default, or if 0 is specified, all available threads are used.

```python
from cogent3 import load_aligned_seqs
from piqtree import jc_distances

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

# Use only 4 threads
distance_matrix = jc_distances(aln, num_threads=4)
```

## See also

- For using the JC distance matrix to construct a rapid neighbour-joining tree, see ["Construct a rapid neighbour-joining tree from a distance matrix"](construct_nj_tree.md).
