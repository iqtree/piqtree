# %% [markdown]
# We evaluate a support for a specific phylogeny using `piqtree_fit`.
# For this simple case, we will assess the support for a specific phylogeny using the GTR model on one alignment.

# %%
import cogent3

from piqtree import download_dataset

aln_path = download_dataset("example.phy.gz", dest_dir="data")
aln = cogent3.load_aligned_seqs(aln_path, moltype="dna", format="phylip")
aln

# %% [markdown]
# We use a tree corresponding to the data in the `example.phy` file.

# %%
tree_path = download_dataset("example.tree.gz", dest_dir="data")
tree = cogent3.load_tree(tree_path)
tree.get_figure().show()

# %% [markdown]
# We now take a look at the help for the `piqtree_fit` app.

# %%

from cogent3 import app_help

app_help("piqtree_fit")

# %% [markdown]
# We fit a GTR model and estimate the nucleotide frequencies by maximum likelihood.

# %%
from cogent3 import get_app

fit_gtr = get_app("piqtree_fit", tree, model="GTR+FO")
fit_gtr

# %% [markdown]
# ## Fit the model

# %%
fit = fit_gtr(aln)

# %% [markdown]
# The `fit` object is a `cogent3` tree (`PhyloNode`) and the maximum likelihood estimated parameters are stored in the `params` attribute.
#
# > **Note**
# > "lnL" is short for log likelihood and is the log likelihood of the model.
# > "mprobs" is short for motif probabilities and are the estimated equilibrium frequencies of the nucleotides.

# %%
fit.params
