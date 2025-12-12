# Find the model of best fit with ModelFinder

IQ-TREE's ModelFinder can be used to automatically find the model of best fit for an alignment using [`model_finder`](../api/model/model_finder.md).
The best scoring model under either **the *Akaike information criterion* (AIC), *corrected Akaike information criterion* (AICc), or the *Bayesian information criterion* (BIC) can be selected.

## Usage

### Basic Usage

Construct a `cogent3` alignment object, then construct a maximum-likelihood tree.

```python
from cogent3 import load_aligned_seqs
from piqtree import model_finder

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

result = model_finder(aln)

best_aic_model = result.best_aic
best_aicc_model = result.best_aicc
best_bic_model = result.best_bic
```

### Specifying the Search Space

We expose the `mset`, `mfreq` and `mrate` parameters from IQ-TREE's ModelFinder to specify the substitution model search space, base frequency search space, and rate heterogeneity search space respectively. They can be specified as a set of strings in either `model_set`, `freq_set` or `rate_set` respectively.

```python
from cogent3 import load_aligned_seqs
from piqtree import model_finder

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

result = model_finder(aln, model_set={"HKY", "TIM"})

best_aic_model = result.best_aic
best_aicc_model = result.best_aicc
best_bic_model = result.best_bic
```

### Reproducible Results

For reproducible results, a random seed may be specified.
> **Caution:** 0 is a specific random seed. None is equivalent to no random seed being specified.

```python
from cogent3 import load_aligned_seqs
from piqtree import model_finder

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

result = model_finder(aln, rand_seed=5)

best_aic_model = result.best_aic
best_aicc_model = result.best_aicc
best_bic_model = result.best_bic
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
from piqtree import model_finder

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

result = model_finder(aln, num_threads=4)

best_aic_model = result.best_aic
best_aicc_model = result.best_aicc
best_bic_model = result.best_bic
```

### Additional options

Additional options in the format that would be passed to the IQ-TREE CLI can be set.
Options in `other_options` already allowed through `model_finder` will be ignored.

```python
from cogent3 import load_aligned_seqs
from piqtree import model_finder

aln = load_aligned_seqs("my_alignment.fasta", moltype="dna")

result = model_finder(aln, model_set={"HKY", "TIM"}, other_options="-mtree")
```

## See also

- For constructing a maximum likelihood tree, see ["Construct a maximum likelihood phylogenetic tree"](construct_ml_tree.md).
- For how to specify a `Model`, see ["Use different kinds of substitution models"](using_substitution_models.md).
