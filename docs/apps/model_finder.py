# %% [markdown]
# We use the `piq_model_finder` app to rank models. This is the python binding to the IQ-TREE ModelFinder tool.

# %%
from cogent3 import app_help, get_app, load_aligned_seqs

from piqtree import download_dataset

aln_path = download_dataset("example.phy.gz", dest_dir="data")
aln = load_aligned_seqs(aln_path, moltype="dna", format="phylip")

# %% [markdown]
# Get help and then apply `piq_model_finder`.

# %%
app_help("piq_model_finder")

# %%
mfinder = get_app("piq_model_finder")
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
# > The process is the same for both the `piq_build_tree` and the `piq_fit_tree` apps.

# %%
fit = get_app("piq_build_tree", selected)
fitted = fit(aln)
fitted
