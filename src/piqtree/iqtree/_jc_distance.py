from collections.abc import Sequence

import numpy as np
from _piqtree import iq_jc_distances
from cogent3.core.alignment import Alignment
from cogent3.evolve.fast_distance import DistanceMatrix

from piqtree.iqtree._decorator import iqtree_func

iq_jc_distances = iqtree_func(iq_jc_distances, hide_files=True)


def _dists_to_distmatrix(
    distances: np.ndarray,
    names: Sequence[str],
) -> DistanceMatrix:
    """Convert numpy representation of distance matrix into cogent3 pairwise distance matrix.

    Parameters
    ----------
    distances : np.ndarray
        Pairwise distances.
    names : Sequence[str]
        Corresponding sequence names.

    Returns
    -------
    DistanceMatrix
        Pairwise distance matrix.

    """
    dist_dict = {}
    for i in range(1, len(distances)):
        for j in range(i):
            dist_dict[(names[i], names[j])] = distances[i, j]
    return DistanceMatrix(dist_dict)


def jc_distances(
    aln: Alignment,
    num_threads: int | None = None,
) -> DistanceMatrix:
    """Compute pairwise JC distances for a given alignment.

    Parameters
    ----------
    aln : Alignment
        The alignment to compute pairwise JC distances for.
    num_threads: int | None, optional
        Number of threads for IQ-TREE to use,
        by default None (uses all available threads).

    Returns
    -------
    DistanceMatrix
        Pairwise JC distance matrix.

    """
    if num_threads is None:
        num_threads = 0

    names = aln.names
    seqs = [str(seq) for seq in aln.iter_seqs(names)]

    distances = np.array(iq_jc_distances(names, seqs, num_threads)).reshape(
        (len(names), len(names)),
    )
    return _dists_to_distmatrix(distances, names)
