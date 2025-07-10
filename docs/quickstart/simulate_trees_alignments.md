# Simulate trees and alignments using msprime and AliSim

A set of coalescent trees and their corresponding alignments can be simulated using `sim_ancestry` from `msprime` and `simulate_alignment` from the `piqtree`'s version of `AliSim`.

## Usage

### Basic Example

The following pipeline simulates a set of trees independently under the coalescent model with a given population size and recombination rate using `msprime`, then rescales the branch lengths of the trees to substitution units, and finally simulates a set of alignments for the rescaled trees using the `piqtree`'s version of `AliSim` with a given substitution model and random seed.

The example below shows the pipeline for generating one pair of tree/alignment for 10 taxa.

```python
from piqtree import simulate_alignment
import cogent3
import msprime

# Sets simulation parameters
num_taxa = 10
seq_length = 2000
recombination_rate = 1e-8
population_size = 10000
num_threads = 5
seed = 1

# Simulates trees under the coalescent model using msprime
ts = msprime.sim_ancestry(
        samples=num_taxa/2,
        sequence_length=seq_length,
        recombination_rate=recombination_rate,
        population_size=population_size,
        random_seed=seed)

for t in ts.trees():
    newick_tree = t.as_newick()
    break
tree = cogent3.make_tree(newick_tree)

# Rescales branch lengths to substitution units
scaling_factor = 1.0 / (2 * population_size)
for node in tree.preorder():
    if node.length is not None:
        node.length *= scaling_factor

print(tree.get_newick(with_distances=True, semicolon=True)) # Prints the rescaled simulated tree

# Simulates alignment from the rescaled simulated tree using AliSim
trees = [tree]
res = simulate_alignment(
        trees = trees,
        subst_model = "JC",
        seed = seed,
        num_threads = num_threads,
        seq_length = seq_length)

print(res[0]) # Prints the alignment
print(res[1]) # Prints logs

```

### Description of Parameters for Tree Simulation

The `sim_ancestry` function in `msprime` can simulate trees under different population genetics models, including the coalescent model with recombination. We explain the parameters used above as well as other useful parameters:
- `samples`: int | dict. Represents the number of sampled individuals in a population.
- `recombination_rate`: float | None. Uniform and non-uniform rates of recombination along the genome.
- `population_size`: int | None. Sets the size of the single constant sized population.
- `ploidy`: int | None. Sets the number of nodes per sample individual, as well as the time scale for continuous time coalescent models. The default value of ploidy is 2, assuming diploid populations.
- `sequence_length`: int | None. Determines the total length of the sequence.
- `model`: str | None. Determines the model under which the ancestral history of the sample is generated. The default model is standard coalescent, but other models such as [Discrete Time Wright-Fisher].(https://tskit.dev/msprime/docs/stable/api.html#msprime.DiscreteTimeWrightFisher) are supported.
- `num_replicates`: int | None. Number of independent simulation replicates (e.g., gene trees).

For further information about the usage of these parameters, see https://tskit.dev/msprime/docs/stable/ancestry.html.

### Description of Parameters for Alignment Simulation

In addition to the input set of trees, `simulate_alignment` allows for specifying several parameters:
- `subst_model`: str | None. Sequence substitution model.
- `num_threads`: int | None. Number of threads (by default None and will be set to 1).
- `partition_info`: list[str] | None. Partition information (by default None and will be set to []).
- `partition_type`: str | None. If provided, partition type must be ‘equal’, ‘proportion’, or ‘unlinked’ (by default None and will be set to "").
- `seq_length`: int | None. The length of sequences (by default None and will be set to 1000).
- `insertion_rate`: float | None. The insertion rate (by default None and will be set to 0.0).
- `deletion_rate`: float | None. The deletion rate (by default None and will be set to 0.0).
- `root_seq`: str | None. The root sequence (by default None and will be set to "").
- `insertion_size_distribution`: str | None. The insertion size distribution (by default None and will be set to "").
- `deletion_size_distribution`: str | None. The deletion size distribution (by default None and will be set to "").

For how to specify a `Model`, see ["Use different kinds of substitution models"](using_substitution_models.md).
