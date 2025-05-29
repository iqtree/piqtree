# %% [markdown]
# We use the `piqtree_mfinder` app to rank models. This is the python binding to the IQ-TREE ModelFinder tool.

# %%
from cogent3 import app_help, get_app, load_aligned_seqs

from piqtree import download_dataset

aln_path = download_dataset("example.phy.gz", dest_dir="data")
# format argument not required after cogent3 2024.11 release
aln = load_aligned_seqs(aln_path, moltype="dna", format="phylip")

# %% [markdown]
# Get help and then apply `piqtree_mfinder`.

# %%
app_help("piqtree_mfinder")

# %%
mfinder = get_app("piqtree_mfinder")
ranked = mfinder(aln)
ranked

# %% [markdown]
# ## Accessing the best model
# The different measures used to select the best model, AIC, AICc, and BIC, are available as attributes of the result object. We'll select AICc as the measure for choosing the best model.

# %%
selected = ranked.best_aicc

# %% [markdown]
# You can inspect the statistics for one of these using the `model_stats` attribute.

# %%
ranked.model_stats[selected]

# %% [markdown]
# ## Using the best model
# You can apply the selected model to a phylogenetic analysis.
# > **Note**
# > The process is the same for both the `piqtree_phylo` and the `piqtree_fit` apps.

# %%
fit = get_app(
    "piqtree_phylo",
    selected.submod_type,
    freq_type=selected.freq_type,
    rate_model=selected.rate_model,
    invariant_sites=selected.invariant_sites,
)
fitted = fit(aln)
fitted
