"""Python wrappers to tree searching functions in the IQ-TREE library."""

from collections.abc import Iterable, Sequence

import cogent3
import cogent3.app.typing as c3_types
import numpy as np
import yaml
from _piqtree import iq_build_tree, iq_consensus_tree, iq_fit_tree, iq_nj_tree
from cogent3 import PhyloNode, make_tree

from piqtree.exceptions import ParseIqTreeError
from piqtree.iqtree._decorator import iqtree_func
from piqtree.model import Model, StandardDnaModel, make_model
from piqtree.util import get_newick

iq_build_tree = iqtree_func(iq_build_tree, hide_files=True)
iq_fit_tree = iqtree_func(iq_fit_tree, hide_files=True)
iq_nj_tree = iqtree_func(iq_nj_tree, hide_files=True)
iq_consensus_tree = iqtree_func(iq_consensus_tree, hide_files=True)

# the order defined in IQ-TREE
RATE_PARS = "A/C", "A/G", "A/T", "C/G", "C/T", "G/T"
MOTIF_PARS = "A", "C", "G", "T"


def _rename_iq_tree(tree: cogent3.PhyloNode, names: Sequence[str]) -> None:
    for tip in tree.tips():
        tip.name = names[int(tip.name)]


def _insert_edge_pars(tree: cogent3.PhyloNode, **kwargs: dict) -> None:
    # inserts the edge parameters into each edge to match the structure of
    # cogent3.PhyloNode
    for node in tree.get_edge_vector():
        # skip the rate parameters when node is the root
        if node.is_root():
            kwargs = {k: v for k, v in kwargs.items() if k == "mprobs"}
            del node.params["edge_pars"]
        node.params.update(kwargs)


def _edge_pars_for_cogent3(tree: cogent3.PhyloNode, model: Model) -> None:
    base_model = model.submod_type.base_model

    rate_pars = tree.params["edge_pars"]["rates"]
    motif_pars = {"mprobs": tree.params["edge_pars"]["mprobs"]}
    # renames parameters to conform to cogent3's naming conventions
    if base_model in {StandardDnaModel.JC, StandardDnaModel.F81}:
        # skip rate_pars since rate parameters are constant in JC and F81
        _insert_edge_pars(
            tree,
            **motif_pars,
        )
        return
    if base_model in {StandardDnaModel.K80, StandardDnaModel.HKY}:
        rate_pars = {"kappa": rate_pars["A/G"]}

    elif base_model is StandardDnaModel.TN:
        rate_pars = {"kappa_r": rate_pars["A/G"], "kappa_y": rate_pars["C/T"]}

    elif base_model is StandardDnaModel.GTR:
        del rate_pars["G/T"]

    # applies global rate parameters to each edge
    _insert_edge_pars(
        tree,
        **rate_pars,
        **motif_pars,
    )


def _parse_nonlie_model(tree: cogent3.PhyloNode, tree_yaml: dict) -> None:
    # parse motif and rate parameters in the tree_yaml for non-Lie DnaModel
    model_fits = tree_yaml.get("ModelDNA", {})

    state_freq_str = model_fits.get("state_freq", "")
    rate_str = model_fits.get("rates", "")

    # parse motif parameters, assign each to a name, and raise an error if not found
    if state_freq_str:
        # converts the strings of motif parameters into dictionary
        state_freq_list = [
            float(value) for value in state_freq_str.replace(" ", "").split(",")
        ]
        tree.params["edge_pars"] = {
            "mprobs": dict(zip(MOTIF_PARS, state_freq_list, strict=True)),
        }
    else:
        msg = "IQ-TREE output malformated, motif parameters not found."
        raise ParseIqTreeError(msg)

    # parse rate parameters, assign each to a name, and raise an error if not found
    if rate_str:
        rate_list = [float(value) for value in rate_str.replace(" ", "").split(",")]
        tree.params["edge_pars"]["rates"] = dict(
            zip(RATE_PARS, rate_list, strict=True),
        )
    else:
        msg = "IQ-TREE output malformated, rate parameters not found."
        raise ParseIqTreeError(msg)


def _parse_lie_model(
    tree: cogent3.PhyloNode,
    tree_yaml: dict,
    lie_model_name: str,
) -> None:
    # parse motif and rate parameters in the tree_yaml for Lie DnaModel
    model_fits = tree_yaml.get(lie_model_name, {})

    # parse motif parameters, assign each to a name, and raise an error if not found
    state_freq_str = model_fits.get("state_freq", "")
    if state_freq_str:
        state_freq_list = [
            float(value) for value in state_freq_str.replace(" ", "").split(",")
        ]
        tree.params[lie_model_name] = {
            "mprobs": dict(zip(MOTIF_PARS, state_freq_list, strict=True)),
        }
    else:
        msg = "IQ-TREE output malformated, motif parameters not found."
        raise ParseIqTreeError(msg)

    # parse rate parameters, skipping LIE_1_1 (aka JC69) since its rate parameter is constant thus absent
    if "model_parameters" in model_fits:
        model_parameters = model_fits["model_parameters"]

        # convert model parameters to a list of floats if they are a string
        if isinstance(model_parameters, str):
            tree.params[lie_model_name]["model_parameters"] = [
                float(value) for value in model_parameters.replace(" ", "").split(",")
            ]
        else:
            # directly use the float
            tree.params[lie_model_name]["model_parameters"] = model_parameters


def _tree_equal(node1: PhyloNode, node2: PhyloNode) -> bool:
    children_group1 = node1.children
    children_group2 = node2.children

    if len(children_group1) != len(children_group2):
        return False

    # recursively check if two PhyloNode have the same name and branch length, so for their children.
    for child1, child2 in zip(children_group1, children_group2, strict=True):
        if not _tree_equal(child1, child2):
            return False

    # handle empty/different internal node names
    if children_group1 == []:
        return node1.name == node2.name and node1.length == node2.length
    return node1.length == node2.length


def _process_tree_yaml(
    tree_yaml: dict,
    names: Sequence[str],
) -> cogent3.PhyloNode:
    newick = tree_yaml["PhyloTree"]["newick"]

    tree = cogent3.make_tree(newick)
    candidates = tree_yaml["CandidateSet"]
    likelihood = None
    for candidate in candidates.values():
        candidate_likelihood, candidate_newick = candidate.split(" ")
        candidate_tree = cogent3.make_tree(candidate_newick)
        if _tree_equal(candidate_tree, tree):
            likelihood = float(candidate_likelihood)
            break
    if likelihood is None:
        msg = "IQ-TREE output malformated, likelihood not found."
        raise ParseIqTreeError(msg)

    tree.params["lnL"] = likelihood

    # parse non-Lie DnaModel parameters
    if "ModelDNA" in tree_yaml:
        _parse_nonlie_model(tree, tree_yaml)

    # parse Lie DnaModel parameters, handling various Lie model names
    elif key := next(
        (key for key in tree_yaml if key.startswith("ModelLieMarkov")),
        None,
    ):
        _parse_lie_model(tree, tree_yaml, key)

    # parse rate model, handling various rate model names
    if key := next((key for key in tree_yaml if key.startswith("Rate")), None):
        tree.params[key] = tree_yaml[key]

    _rename_iq_tree(tree, names)

    tree.name_unnamed_nodes()

    return tree


def build_tree(
    aln: c3_types.AlignedSeqsType,
    model: Model | str,
    rand_seed: int | None = None,
    bootstrap_replicates: int | None = None,
    num_threads: int | None = None,
) -> cogent3.PhyloNode:
    """Reconstruct a phylogenetic tree.

    Given a sequence alignment, uses IQ-TREE to reconstruct a phylogenetic tree.

    Parameters
    ----------
    aln : c3_types.AlignedSeqsType
        The sequence alignment.
    model : Model | str
        The substitution model with base frequencies and rate heterogeneity.
    rand_seed : int | None, optional
        The random seed - 0 or None means no seed, by default None.
    bootstrap_replicates : int, optional
        The number of bootstrap replicates to perform, by default None.
        If 0 is provided, then no bootstrapping is performed.
        At least 1000 is required to perform bootstrapping.
    num_threads: int | None, optional
        Number of threads for IQ-TREE to use, by default None (single-threaded).
        If 0 is specified, IQ-TREE attempts to find the optimal number of threads.

    Returns
    -------
    cogent3.PhyloNode
        The IQ-TREE maximum likelihood tree from the given alignment.

    """
    if isinstance(model, str):
        model = make_model(model)

    if rand_seed is None:
        rand_seed = 0  # The default rand_seed in IQ-TREE

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
        ),
    )
    tree = _process_tree_yaml(yaml_result, names)

    # for non-Lie models, populate parameters to each branch and
    # rename them to mimic cogent3.PhyloNode
    if "edge_pars" in tree.params:
        _edge_pars_for_cogent3(tree, model)
    return tree


def fit_tree(
    aln: c3_types.AlignedSeqsType,
    tree: cogent3.PhyloNode,
    model: Model | str,
    num_threads: int | None = None,
    *,
    bl_fixed: bool = False,
) -> cogent3.PhyloNode:
    """Fit branch lengths and likelihood for a tree.

    Given a sequence alignment and a fixed topology,
    uses IQ-TREE to fit branch lengths to the tree.

    Parameters
    ----------
    aln : c3_types.AlignedSeqsType
        The sequence alignment.
    tree : cogent3.PhyloNode
        The topology to fit branch lengths to.
    model : Model | str
        The substitution model with base frequencies and rate heterogeneity.
    num_threads: int | None, optional
        Number of threads for IQ-TREE to use, by default None (single-threaded).
        If 0 is specified, IQ-TREE attempts to find the optimal number of threads.
    bl_fixed: bool, optional
        If True, evaluates likelihood using the provided branch lengths on the tree.
        Branch lengths will be treated as constant in this case. Otherwise if False,
        branch lengths are fitted to the tree whether provided or not. By default False.

    Returns
    -------
    cogent3.PhyloNode
        A phylogenetic tree with same given topology fitted with branch lengths.

    """
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
        ),
    )
    tree = _process_tree_yaml(yaml_result, names)

    # for non-Lie models, populate parameters to each branch and
    # rename them to mimic cogent3.PhyloNode
    if "edge_pars" in tree.params:
        _edge_pars_for_cogent3(tree, model)
    return tree


def nj_tree(
    pairwise_distances: c3_types.PairwiseDistanceType,
    *,
    allow_negative: bool = False,
) -> cogent3.PhyloNode:
    """Construct a neighbour joining tree from a pairwise distance matrix.

    Parameters
    ----------
    pairwise_distances : c3_types.PairwiseDistanceType
        Pairwise distances to construct neighbour joining tree from.
    allow_negative : bool, optional
        Whether to allow negative branch lengths in the output.
        Coerces to 0 if not allowed, by default False.

    Returns
    -------
    cogent3.PhyloNode
        The neigbour joining tree.

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
            node.length = max(node.length, 0)

    return tree


def _all_same_taxa_set(trees: Iterable[cogent3.PhyloNode]) -> bool:
    tree_it = iter(trees)
    try:
        taxa_set = set(next(tree_it).get_tip_names())
    except StopIteration:
        return True

    return all(taxa_set == set(tree.get_tip_names()) for tree in tree_it)


def consensus_tree(
    trees: Iterable[cogent3.PhyloNode],
    *,
    min_support: float = 0.5,
) -> cogent3.PhyloNode:
    """Build a consensus tree, defaults to majority-rule consensus tree.

    The min_support parameter represents the proportion of trees a clade
    must appear in to be in the resulting consensus tree.

    If min_support is 1.0, computes the strict consensus tree.
    If min_support is 0.0, computes the extended majority-rule consensus tree.

    Parameters
    ----------
    trees : Iterable[cogent3.PhyloNode]
        The trees to form a consensus tree from.
    min_support : float, optional
        The minimum support for a clade to appear
        in the consensus tree, by default 0.5.

    Returns
    -------
    cogent3.PhyloNode
        The constructed consensus tree.

    """
    if not 0 <= min_support <= 1:
        msg = f"Only min support values in the range 0 <= value < 1 are supported, got {min_support}"
        raise ValueError(msg)

    if not _all_same_taxa_set(trees):
        msg = "Trees must be on same taxa set."
        raise ValueError(msg)

    newick_trees = [get_newick(tree) for tree in trees]

    return make_tree(iq_consensus_tree(newick_trees, min_support))
