# %% [markdown]
# ## Building workflows using `piqtree` apps
# > **WARNING**
# > This page is under construction!
#
# We can combine `piqtree` apps with other cogent3 apps to develop a pipeline. There are multiple concepts involved here, particularly data stores, composed apps, parallel execution, log files etc... See the `cogent3` [app documentation](https://cogent3.org/doc/app/index.html) for more details.
#
# To develop a pipeline efficiently we only need a subset of the sequences in an alignment. We will use the [diverse-seq](https://pypi.org/project/diverse-seq/) plugin for that purpose. This allows selecting a specified subset of sequences that capture the diversity in an alignment.
#
# But first, we need the data.
# %%
from piqtree import download_dataset

alns_path = download_dataset("mammal-orths.zip", dest_dir="data", inflate_zip=False)

# %% [markdown]
# We open this directory as a `cogent3` data store.

# %%
from cogent3 import open_data_store

dstore = open_data_store(alns_path, suffix="fa")
dstore.describe

# %% [markdown]
# We need to create some apps to: load data, a divergent sequence selector app, drop alignment columns containing non-canonical nucleotides (so gaps and N's), select alignments with a minimym number of aligned columns, a "data store" to write results to and a writer. These will then be combined into a single composed app which will be applied to all the alignments in the data store.

# %%
import pathlib
from collections import Counter

from cogent3 import get_app

outpath = pathlib.Path("data/delme.sqlitedb")
outpath.unlink(missing_ok=True)

loader = get_app("load_aligned", format="fasta", moltype="dna")
divergent = get_app("dvs_nmost", n=10, k=6)
just_nucs = get_app("omit_degenerates")  # has to go after the divergent selector
min_length = get_app("min_length", length=600)
best_model = get_app("piq_model_finder")
app = loader + divergent + just_nucs + min_length + best_model
model_counts = Counter(
    str(result.best_aic)
    for result in app.as_completed(dstore, show_progress=True, parallel=True)
    if result
)
model_counts
