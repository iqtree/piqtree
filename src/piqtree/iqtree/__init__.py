"""Functions for calling IQ-TREE as a library."""

from ._alignment import simulate_alignment
from ._fit_tree_hessian import HessianFit, fit_tree_hessian
from ._jc_distance import jc_distances
from ._model_finder import ModelFinderResult, ModelResultValue, model_finder
from ._random_tree import TreeGenMode, random_tree
from ._robinson_foulds import robinson_foulds
from ._tree import build_tree, consensus_tree, fit_tree, nj_tree

__all__ = [
    "HessianFit",
    "ModelFinderResult",
    "ModelResultValue",
    "TreeGenMode",
    "build_tree",
    "consensus_tree",
    "fit_tree",
    "fit_tree_hessian",
    "jc_distances",
    "model_finder",
    "nj_tree",
    "random_tree",
    "robinson_foulds",
    "simulate_alignment",
]
