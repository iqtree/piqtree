"""Python wrappers to RF distances in the IQ-TREE library."""

from collections.abc import Sequence

import numpy as np
from _piqtree import iq_robinson_fould
from cogent3.core.tree import PhyloNode

from piqtree.iqtree._decorator import iqtree_func
from piqtree.util import get_newick

iq_robinson_fould = iqtree_func(iq_robinson_fould)


def robinson_foulds(trees: Sequence[PhyloNode]) -> np.ndarray:
    """Pairwise Robinson-Foulds distance between a sequence of trees.

    For the given collection of trees, returns a numpy array containing
    the pairwise distances between the trees.

    Parameters
    ----------
    trees : Sequence[PhyloNode]
        The sequence of trees to calculate the pairwise Robinson-Foulds
        distances of.

    Returns
    -------
    np.ndarray
        Pairwise Robinson-Foulds distances.

    """
    pairwise_distances = np.zeros((len(trees), len(trees)))
    for i in range(1, len(trees)):
        for j in range(i):
            rf = iq_robinson_fould(get_newick(trees[i]), get_newick(trees[j]))
            pairwise_distances[i, j] = rf
            pairwise_distances[j, i] = rf
    return pairwise_distances
