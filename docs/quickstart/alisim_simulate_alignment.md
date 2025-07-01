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

## See also

- For how to specify a `Model`, see ["Use different kinds of substitution models"](using_substitution_models.md).
