"""piqtree - access the power of IQ-TREE within Python."""


def _add_dll_path() -> None:
    import os
    from pathlib import Path

    dll_folder = Path(__file__).parent / "_libiqtree"
    if "add_dll_directory" in dir(os) and dll_folder.is_dir():
        os.add_dll_directory(dll_folder)  # type: ignore[attr-defined]


_add_dll_path()
del _add_dll_path

from piqtree._data import dataset_names, download_dataset
from piqtree._piqtree import __iqtree_version__
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
    simulate_alignment,
)
from piqtree.model import (
    Model,
    available_freq_type,
    available_models,
    available_rate_type,
    make_model,
)

__version__ = "0.8.4"


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
    "simulate_alignment",
]
