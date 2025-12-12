"""Python wrappers to tree searching functions in the IQ-TREE library."""

from collections.abc import Iterable, Sequence
from typing import Any, cast

import numpy as np
import yaml
from _piqtree import iq_build_tree, iq_consensus_tree, iq_fit_tree, iq_nj_tree
from cogent3 import make_tree
from cogent3.core.alignment import Alignment
from cogent3.core.tree import PhyloNode
from cogent3.evolve.fast_distance import DistanceMatrix

from piqtree.exceptions import ParseIqTreeError
from piqtree.iqtree._decorator import iqtree_func
from piqtree.iqtree._parse_tree_parameters import parse_model_parameters
from piqtree.model import Model, make_model
from piqtree.util import get_newick, process_rand_seed_nonzero, validate_other_options

iq_build_tree = iqtree_func(iq_build_tree, hide_files=True)
iq_fit_tree = iqtree_func(iq_fit_tree, hide_files=True)
iq_nj_tree = iqtree_func(iq_nj_tree, hide_files=True)
iq_consensus_tree = iqtree_func(iq_consensus_tree, hide_files=True)


def _rename_iq_tree(tree: PhyloNode, names: Sequence[str]) -> None:
    for tip in tree.tips():
        tip.name = names[int(tip.name)]


def _tree_equal(node1: PhyloNode, node2: PhyloNode) -> bool:
    children_group1 = node1.children
    children_group2 = node2.children

    if len(children_group1) != len(children_group2):
        return False

    # recursively check if two PhyloNodes have the same name and branch length, and do the same for their children.
    for child1, child2 in zip(children_group1, children_group2, strict=True):
        if not _tree_equal(child1, child2):
            return False

    # handle empty/different internal node names
    if children_group1 == []:
        return node1.name == node2.name and node1.length == node2.length
    return node1.length == node2.length


def _process_tree_yaml(
    tree_yaml: dict[str, Any],
    names: Sequence[str],
    model: Model,
) -> PhyloNode:
    newick = tree_yaml["PhyloTree"]["newick"]

    tree = make_tree(newick)
    candidates = tree_yaml["CandidateSet"]
    likelihood = None
    for candidate in candidates.values():
        candidate_likelihood, candidate_newick = candidate.split(" ")
        candidate_tree = make_tree(candidate_newick)
        if _tree_equal(candidate_tree, tree):
            likelihood = float(candidate_likelihood)
            break
    if likelihood is None:
        msg = "IQ-TREE output is malformed, likelihood not found."
        raise ParseIqTreeError(msg)

    tree.params["lnL"] = likelihood

    parse_model_parameters(tree, tree_yaml, model)

    # parse rate model, handling various rate model names
    if key := next((key for key in tree_yaml if key.startswith("Rate")), None):
        tree.params[key] = tree_yaml[key]

    _rename_iq_tree(tree, names)

    tree.name_unnamed_nodes()
    tree.params["model"] = str(model)

    return tree


# Bad options
INVALID_BUILD_TREE_PARAMS = [
    "-s",  # aln file
    "-m",  # model selection
    "-seed",  # seed
    "-bb",  # bootstrap replicates
    "-nt",  # threads
    "-ntmax",  # threads
]


def build_tree(
    aln: Alignment,
    model: Model | str,
    rand_seed: int | None = None,
    bootstrap_replicates: int | None = None,
    num_threads: int | None = None,
    other_options: str = "",
) -> PhyloNode:
    """Reconstruct a phylogenetic tree.

    Given a sequence alignment, uses IQ-TREE to reconstruct a phylogenetic tree.

    Parameters
    ----------
    aln : Alignment
        The sequence alignment.
    model : Model | str
        The substitution model with base frequencies and rate heterogeneity.
    rand_seed : int | None, optional
        The random seed - None means no seed is used, by default None.
    bootstrap_replicates : int | None, optional
        The number of bootstrap replicates to perform, by default None.
        If 0 is provided, then no bootstrapping is performed.
        At least 1000 is required to perform bootstrapping.
    num_threads: int | None, optional
        Number of threads for IQ-TREE to use, by default None (single-threaded).
        If 0 is specified, IQ-TREE attempts to find the optimal number of threads.
    other_options: str, optional
        Additional command line options for IQ-TREE.

    Returns
    -------
    PhyloNode
        The IQ-TREE maximum likelihood tree from the given alignment.

    """
    validate_other_options(other_options, INVALID_BUILD_TREE_PARAMS)

    if isinstance(model, str):
        model = make_model(model)

    rand_seed = process_rand_seed_nonzero(rand_seed)

    if bootstrap_replicates is None:
        bootstrap_replicates = 0

    if num_threads is None:
        num_threads = 1

    names = aln.names
    seqs = [str(seq) for seq in aln.iter_seqs(names)]

    yaml_result = yaml.safe_load(
        iq_build_tree(
            names,
            seqs,
            str(model),
            rand_seed,
            bootstrap_replicates,
            num_threads,
            other_options,
        ),
    )
    return _process_tree_yaml(yaml_result, names, model)


INVALID_FIT_TREE_PARAMS = [
    "-s",  # aln file
    "-t",  # tree specification
    "-te",  # tree specification
    "-m",  # model selection
    "-nt",  # threads
    "-ntmax",  # threads
    "-blfix",  # whether to fix current branch lengths
]


def fit_tree(
    aln: Alignment,
    tree: PhyloNode,
    model: Model | str,
    num_threads: int | None = None,
    other_options: str = "",
    *,
    bl_fixed: bool = False,
) -> PhyloNode:
    """Fit branch lengths and likelihood for a tree.

    Given a sequence alignment and a fixed topology,
    uses IQ-TREE to fit branch lengths to the tree.

    Parameters
    ----------
    aln : Alignment
        The sequence alignment.
    tree : PhyloNode
        The topology to fit branch lengths to.
    model : Model | str
        The substitution model with base frequencies and rate heterogeneity.
    num_threads: int | None, optional
        Number of threads for IQ-TREE to use, by default None (single-threaded).
        If 0 is specified, IQ-TREE attempts to find the optimal number of threads.
    bl_fixed: bool, optional
        If True, evaluates likelihood using the provided branch lengths on the tree.
        Branch lengths will be treated as constant in this case, with any unspecified
        branch lengths still being optimised. Otherwise if False, branch lengths are
        fitted to the tree whether provided or not. By default False.
    other_options: str, optional
        Additional command line options for IQ-TREE.

    Returns
    -------
    PhyloNode
        A phylogenetic tree with the same given topology fitted with branch lengths.

    """
    validate_other_options(other_options, INVALID_FIT_TREE_PARAMS)

    if isinstance(model, str):
        model = make_model(model)

    if num_threads is None:
        num_threads = 1

    names = aln.names
    seqs = [str(seq) for seq in aln.iter_seqs(names)]
    newick = get_newick(tree)

    yaml_result = yaml.safe_load(
        iq_fit_tree(
            names,
            seqs,
            str(model),
            newick,
            bl_fixed,
            0,
            num_threads,
            other_options,
        ),
    )
    return _process_tree_yaml(yaml_result, names, model)


def nj_tree(
    pairwise_distances: DistanceMatrix,
    *,
    allow_negative: bool = False,
) -> PhyloNode:
    """Construct a neighbour joining tree from a pairwise distance matrix.

    Parameters
    ----------
    pairwise_distances : DistanceMatrix
        Pairwise distances to construct neighbour joining tree from.
    allow_negative : bool, optional
        Whether to allow negative branch lengths in the output.
        Coerces to 0 if not allowed, by default False.

    Returns
    -------
    PhyloNode
        The neighbour joining tree.

    See Also
    --------
    jc_distances : construction of pairwise JC distance matrix from alignment.

    """
    if np.isnan(pairwise_distances.array).any():
        msg = "The pairwise distance matrix cannot contain NaN values."
        raise ValueError(msg)

    newick_tree = iq_nj_tree(
        pairwise_distances.keys(),
        np.array(pairwise_distances).flatten(),
    )

    tree = make_tree(newick_tree)

    if not allow_negative:
        for node in tree.preorder(include_self=False):
            node.length = max(cast("float", node.length), 0)

    return tree


def _all_same_taxa_set(trees: Iterable[PhyloNode]) -> bool:
    tree_it = iter(trees)
    try:
        taxa_set = set(next(tree_it).get_tip_names())
    except StopIteration:
        return True

    return all(taxa_set == set(tree.get_tip_names()) for tree in tree_it)


def consensus_tree(
    trees: Iterable[PhyloNode],
    *,
    min_support: float = 0.5,
) -> PhyloNode:
    """Build a consensus tree, defaults to majority-rule consensus tree.

    The min_support parameter represents the proportion of trees a clade
    must appear in to be in the resulting consensus tree.

    If min_support is 1.0, computes the strict consensus tree.
    If min_support is 0.0, computes the extended majority-rule consensus tree.

    Parameters
    ----------
    trees : Iterable[PhyloNode]
        The trees to form a consensus tree from.
    min_support : float, optional
        The minimum support for a clade to appear
        in the consensus tree, by default 0.5.

    Returns
    -------
    PhyloNode
        The constructed consensus tree.

    """
    if not 0 <= min_support <= 1:
        msg = f"Only min support values in the range 0 <= value <= 1 are supported, got {min_support}"
        raise ValueError(msg)

    if not _all_same_taxa_set(trees):
        msg = "Trees must be on same taxa set."
        raise ValueError(msg)

    newick_trees = [get_newick(tree) for tree in trees]

    return make_tree(iq_consensus_tree(newick_trees, min_support))
