from collections.abc import Sequence

import numpy as np
from cogent3.core.alignment import Alignment
from cogent3.evolve.fast_distance import DistanceMatrix

from piqtree._piqtree import iq_jc_distances
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
    order = sorted(range(len(names)), key=names.__getitem__)
    lower = np.tril(distances, -1)
    symmetric = (lower + lower.T)[np.ix_(order, order)]
    return DistanceMatrix.from_array_names(symmetric, [names[i] for i in order])


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

    distances = np.asarray(
        iq_jc_distances(names=names, seqs=seqs, num_threads=num_threads),
    ).reshape(
        (len(names), len(names)),
    )
    return _dists_to_distmatrix(distances, names)
