"""piqtree - access the power of IQ-TREE within Python."""

from _piqtree import __iqtree_version__

from piqtree._data import dataset_names, download_dataset
from piqtree.iqtree import (
    ModelFinderResult,
    TreeGenMode,
    build_tree,
    fit_tree,
    jc_distances,
    model_finder,
    nj_tree,
    random_trees,
    robinson_foulds,
)
from piqtree.model import (
    Model,
    available_freq_type,
    available_models,
    available_rate_type,
    make_model,
)

__version__ = "0.5.0"

__all__ = [
    "Model",
    "ModelFinderResult",
    "TreeGenMode",
    "__iqtree_version__",
    "available_freq_type",
    "available_models",
    "available_rate_type",
    "build_tree",
    "dataset_names",
    "download_dataset",
    "fit_tree",
    "jc_distances",
    "make_model",
    "model_finder",
    "nj_tree",
    "random_trees",
    "robinson_foulds",
]
