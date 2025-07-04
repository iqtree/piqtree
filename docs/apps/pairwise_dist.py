# %% [markdown]
# We estimate pairwise distances via the Jukes-Cantor model with the `piq_jc_distances` app.

# %%
import cogent3

from piqtree import download_dataset

aln_path = download_dataset("example.phy.gz", dest_dir="data")
aln = cogent3.load_aligned_seqs(aln_path, moltype="dna", format="phylip")

# %% [markdown]
# We get help on the `piqtree_jc_dist` app.

# %%
cogent3.app_help("piq_jc_distances")

# %% [markdown]
# Make an app and apply it to the alignment.

# %%
jc_dists = cogent3.get_app("piq_jc_distances")
dists = jc_dists(aln)
dists
