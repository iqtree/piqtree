"""Fit a tree and return the log-likelihood gradient and Hessian matrix."""

from dataclasses import dataclass
from typing import cast

import numpy as np
from _piqtree import iq_fit_tree_hessian
from cogent3 import make_tree
from cogent3.core.alignment import Alignment
from cogent3.core.tree import PhyloNode

from piqtree.iqtree._decorator import iqtree_func
from piqtree.iqtree._tree import _rename_iq_tree
from piqtree.model import Model, make_model
from piqtree.util import get_newick, process_rand_seed_nonzero, validate_other_options

iq_fit_tree_hessian = iqtree_func(iq_fit_tree_hessian, hide_files=True)


@dataclass(frozen=True)
class HessianFit:
    """Result of :func:`fit_tree_hessian`.

    Attributes
    ----------
    tree : PhyloNode
        The (rooted) fitted tree. The branch order returned by IQ-TREE matches
        ``branch_lengths`` element-for-element.
    branch_lengths : numpy.ndarray
        1-D float64 array of length ``n_branches``.
    gradient : numpy.ndarray
        1-D float64 array of first derivatives of the log-likelihood with
        respect to each branch length, in the same order as ``branch_lengths``.
    hessian : numpy.ndarray
        2-D float64 array of shape ``(n_branches, n_branches)`` — the Hessian
        of the log-likelihood with respect to branch lengths.
    log_likelihood : float
        Log-likelihood at the fitted point.
    """

    tree: PhyloNode
    branch_lengths: np.ndarray
    gradient: np.ndarray
    hessian: np.ndarray
    log_likelihood: float


# Options forbidden because they'd conflict with what this function controls.
_INVALID_FIT_TREE_HESSIAN_PARAMS = [
    "-s",  # aln file
    "-t",  # tree spec
    "-te",  # tree spec
    "-m",  # model selection
    "-nt",  # threads
    "-ntmax",  # threads
    "-blfix",  # fix branch lengths
    "--dating",  # we always pass --dating mcmctree internally
]


def fit_tree_hessian(
    aln: Alignment,
    tree: PhyloNode,
    model: Model | str,
    rand_seed: int | None = None,
    num_threads: int | None = None,
    other_options: str = "",
    *,
    bl_fixed: bool = False,
) -> HessianFit:
    """Fit a tree and return its log-likelihood gradient and Hessian matrix.

    Given a sequence alignment and a topology, IQ-TREE fits branch lengths and
    model parameters, then computes the first derivatives (gradient) and the
    full Hessian matrix of the log-likelihood with respect to the branch
    lengths at the fitted point. This is the same information IQ-TREE would
    normally emit to a ``.mcmctree.hessian`` file via ``--dating mcmctree``,
    but here it is captured in memory without disk I/O.

    Parameters
    ----------
    aln : Alignment
        The sequence alignment.
    tree : PhyloNode
        The topology to fit branch lengths to.
    model : Model | str
        The substitution model with base frequencies and rate heterogeneity.
        Partitioned (super-tree) models are not supported.
    rand_seed : int | None, optional
        The random seed. None (default) means a fresh seed is generated;
        passing the same seed across runs makes the analysis deterministic.
    num_threads : int | None, optional
        Number of threads for IQ-TREE to use, by default None (single-threaded).
        If 0 is specified, IQ-TREE attempts to find the optimal number of
        threads.
    bl_fixed : bool, optional
        If True, evaluates likelihood and derivatives using the provided
        branch lengths on the tree. Otherwise branch lengths are first
        optimised. By default False.
    other_options : str, optional
        Additional command-line options for IQ-TREE.

    Returns
    -------
    HessianFit
        The fitted tree together with branch lengths, gradient and Hessian.

    Notes
    -----
    The branch ordering used by IQ-TREE's Hessian output is the rooted-tree
    ordering described in the MCMCTree documentation. ``branch_lengths``,
    ``gradient`` and the rows/columns of ``hessian`` all share that ordering;
    the corresponding tree topology is returned in ``tree``.
    """
    validate_other_options(other_options, _INVALID_FIT_TREE_HESSIAN_PARAMS)

    if isinstance(model, str):
        model = make_model(model)

    if num_threads is None:
        num_threads = 1

    names = aln.names
    seqs = [str(seq) for seq in aln.iter_seqs(names)]
    newick = get_newick(tree)

    rand_seed = process_rand_seed_nonzero(rand_seed)

    result = cast(
        "dict[str, object]",
        iq_fit_tree_hessian(
            names,
            seqs,
            str(model),
            newick,
            bl_fixed,
            rand_seed,
            num_threads,
            other_options,
        ),
    )

    fitted_tree = make_tree(cast("str", result["tree"]))
    _rename_iq_tree(fitted_tree, names)
    fitted_tree.name_unnamed_nodes()
    fitted_tree.params["model"] = str(model)
    log_likelihood = float(cast("float", result["log_likelihood"]))
    fitted_tree.params["lnL"] = log_likelihood

    return HessianFit(
        tree=fitted_tree,
        branch_lengths=cast("np.ndarray", result["branch_lengths"]),
        gradient=cast("np.ndarray", result["gradient"]),
        hessian=cast("np.ndarray", result["hessian"]),
        log_likelihood=log_likelihood,
    )
