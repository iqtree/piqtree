# %% [markdown]
# The Neighbour-Joining method uses genetic distances to build a phylogenetic tree. `piqtree` provides only `piq_jc_distances` for this. `cogent3` includes many more methods. The  results of either can be used to build a tree. For divergent sequences we will use Lake's paralinear measure as it accomodates divergent sequence compositions.

# %%
import cogent3

from piqtree import download_dataset

aln_path = download_dataset("example.phy.gz", dest_dir="data")
aln = cogent3.load_aligned_seqs(aln_path, moltype="dna", format_name="phylip")

# %% [markdown]
# ## Getting a paralinear distance matrix
# This can be obtained directly from the alignment object itself.

# %%
dists = aln.distance_matrix(calc="paralinear")
dists

# %% [markdown]
# Get help on the `piq_nj_tree` app.

# %%
cogent3.app_help("piq_nj_tree")

# %% [markdown]
# Make an app and apply it to the distance matrix.

# %%
nj = cogent3.get_app("piq_nj_tree")
tree = nj(dists)

# %% [markdown]
# > **Warning**
# > Branch lengths can be negative in the piqtree NJ tree. This manifests as branches going backwards!

# %%
tree.get_figure().show()

# %% [markdown]
# > **Note**
# > To write the tree to a file, use the `write()` method.

# %% [markdown]
# ## Combining the piqtree dist and nj apps
# We can combine the `piq_jc_distances` and `piq_nj_tree` apps to build a tree from an alignment in one step.

# %%
jc = cogent3.get_app("piq_jc_distances")
nj = cogent3.get_app("piq_nj_tree")
app = jc + nj
tree = app(aln)
tree.get_figure().show()
