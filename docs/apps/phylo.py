# %% [markdown]
# We use the `piqtree_phylo` app to build phylogenies.
# For this simple case, we will build a single phylogeny using the GTR model on one alignment.
# We have a utility script for this documentation that provides the local path to that data. We will then load that data and, as it contains quite a few sequences, we will use a subset of it. We use methods on the cogent3 object to do that.

# %%
import cogent3

from piqtree import download_dataset

aln_path = download_dataset("example.phy.gz", dest_dir="data")
aln = cogent3.load_aligned_seqs(aln_path, moltype="dna", format="phylip")
aln


# %% [markdown]
# We now take a look at the help for the `piqtree_phylo` app.

# %%

from cogent3 import app_help

app_help("piqtree_phylo")

# %% [markdown]
# We build an app for estimating phylogenies with a GTR model

# %%
from cogent3 import get_app

phylo_gtr = get_app("piqtree_phylo", model="GTR", bootstrap_reps=1000)
phylo_gtr

# %% [markdown]
# ## Run the phylogeny estimation and display branches with support < 90%.

# %%
tree = phylo_gtr(aln)
dnd = tree.get_figure(show_support=True, threshold=90)
dnd.show(width=600, height=600)
