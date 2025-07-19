"""piqtree - access the power of IQ-TREE within Python."""


def _add_dll_path() -> None:
    import os

    if "add_dll_directory" in dir(os):
        dll_folder = os.path.join(os.path.dirname(__file__), "_libiqtree")
        os.add_dll_directory(dll_folder)  # type: ignore[attr-defined]


_add_dll_path()
del _add_dll_path


from _piqtree import __iqtree_version__

from piqtree._data import dataset_names, download_dataset
from piqtree.iqtree import (
    ModelFinderResult,
    TreeGenMode,
    build_tree,
    consensus_tree,
    fit_tree,
    jc_distances,
    model_finder,
    nj_tree,
    random_tree,
    robinson_foulds,
)
from piqtree.model import (
    Model,
    available_freq_type,
    available_models,
    available_rate_type,
    make_model,
)

__version__ = "0.6.1"

__all__ = [
    "Model",
    "ModelFinderResult",
    "TreeGenMode",
    "__iqtree_version__",
    "available_freq_type",
    "available_models",
    "available_rate_type",
    "build_tree",
    "consensus_tree",
    "dataset_names",
    "download_dataset",
    "fit_tree",
    "jc_distances",
    "make_model",
    "model_finder",
    "nj_tree",
    "random_tree",
    "robinson_foulds",
]
