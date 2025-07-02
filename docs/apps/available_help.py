# %% [markdown]
# To see the apps that `piqtree` makes available, use the `cogent3` function `available_apps()` as follows.

# %%
from cogent3 import available_apps

available_apps("piqtree")

# %% [markdown]
# For a particular app, use the `cogent3` function `app_help()` to find out what the options are.

# %%
from cogent3 import app_help

app_help("piq_build_tree")

# %%
from cogent3 import app_help

app_help("piq_fit_tree")
