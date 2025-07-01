# Execute AliSim simulation and output the alignment

An alignment can be simulated using `simulate_alignment`, given a list of trees and a model

## Usage

### Basic Usage

Construct `cogent3` tree(s), choose the model and random seed, then simulate an alignment using AliSim.

```python
from piqtree import simulate_alignment
import cogent3

tree = cogent3.make_tree("((A:0.1,B:0.2):0.1,(C:0.3,D:0.4):0.2,E:0.5);")
trees = [tree]

res = simulate_alignment(trees, "JC", 1)

print(res[0]) # Prints the alignment
print(res[1]) # Prints logs
```

### Multithreading

To speed up computation, the number of threads to be used may be specified.
By default, the computation is done on a single thread. If 0 is specified,
then IQ-TREE attempts to determine the optimal number of threads.

```python
from piqtree import simulate_alignment
import cogent3

tree = cogent3.make_tree("((A:0.1,B:0.2):0.1,(C:0.3,D:0.4):0.2,E:0.5);")
trees = [tree]

res = simulate_alignment(trees, "JC", 1, num_threads = 4)

print(res[0]) # Prints the alignment
print(res[1]) # Prints logs
```

### Other (Optional) Input Parameters

Apart from the parameters above, simulate_alignment also allows specifying several other parameters described below:
- `partition_info`: partition_info: list[str] | None. Partition information (by default None and will be set to []).
- `partition_type`: str | None. If provided, partition type must be ‘equal’, ‘proportion’, or ‘unlinked’ (by default None and will be set to "").
- `seq_length`: int | None. The length of sequences (by default None and will be set to 1000).
- `insertion_rate`: float | None. The insertion rate (by default None and will be set to 0.0).
- `deletion_rate`: float | None. The deletion rate (by default None and will be set to 0.0).
- `root_seq`: str | None. The root sequence (by default None and will be set to "").
- `insertion_size_distribution`: str | None. The insertion size distribution (by default None and will be set to "").
- `deletion_size_distribution`: str | None. The deletion size distribution (by default None and will be set to "").

## See also

- For how to specify a `Model`, see ["Use different kinds of substitution models"](using_substitution_models.md).
