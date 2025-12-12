from typing import cast

import numpy as np
import pytest
from cogent3 import get_app, make_tree
from cogent3.app.result import model_result
from cogent3.core.alignment import Alignment
from cogent3.core.tree import PhyloNode

import piqtree
from piqtree.exceptions import IqTreeError
from piqtree.model import Model, StandardDnaModel


def check_model_name(got: PhyloNode, expected: str) -> None:
    got_model = cast("str", got.params.get("model"))
    assert got_model == expected


def check_likelihood(got: PhyloNode, expected: model_result) -> None:
    assert got.params["lnL"] == pytest.approx(expected.lnL)


def check_motif_probs(got: PhyloNode, expected: PhyloNode) -> None:
    expected_mprobs = expected.params["mprobs"]
    got_mprobs = got.params["mprobs"]

    expected_keys = set(expected_mprobs.keys())
    got_keys = set(got_mprobs.keys())

    # Check that the base characters are the same
    assert expected_keys == got_keys

    # Check that the probs are the same
    expected_values = [expected_mprobs[key] for key in expected_keys]
    got_values = [got_mprobs[key] for key in expected_keys]
    assert all(
        got == pytest.approx(exp)
        for got, exp in zip(got_values, expected_values, strict=True)
    )


def check_rate_parameters(got: PhyloNode, expected: PhyloNode) -> None:
    # Collect all rate parameters in got and expected
    exclude = {"length", "ENS", "paralinear", "mprobs"}
    expected_keys = {
        k for k in expected.get_edge_vector()[0].params if k not in exclude
    }
    got_keys = {k for k in got.get_edge_vector()[0].params if k not in exclude}

    # Check that the keys of rate are the same
    assert expected_keys == got_keys

    # Check that the values of rate are the same
    expected_values = [expected[0].params[key] for key in expected_keys]
    got_values = [got[0].params[key] for key in expected_keys]

    assert all(
        got == pytest.approx(exp, rel=1e-2)
        for got, exp in zip(got_values, expected_values, strict=True)
    )


def check_branch_lengths(got: PhyloNode, expected: PhyloNode) -> None:
    got_dists = got.tip_to_tip_distances()
    expected_dists = expected.tip_to_tip_distances()
    # make sure the distance matrices have the same name order
    # so we can just compare entire numpy arrays
    expected_dists = expected_dists.take_dists(got_dists.names)
    # Check that the keys of branch lengths are the same
    assert set(got_dists.names) == set(expected_dists.names)

    # Check that the branch lengths are the same
    np.testing.assert_allclose(got_dists.array, expected_dists.array, atol=1e-4)


@pytest.mark.parametrize(
    ("iq_model", "c3_model"),
    [
        (StandardDnaModel.JC, "JC69"),
        (StandardDnaModel.K80, "K80"),
        (StandardDnaModel.GTR, "GTR"),
        (StandardDnaModel.TN, "TN93"),
        (StandardDnaModel.HKY, "HKY85"),
        (StandardDnaModel.F81, "F81"),
    ],
)
def test_fit_tree(
    three_otu: Alignment,
    iq_model: StandardDnaModel,
    c3_model: str,
) -> None:
    tree_topology = make_tree(tip_names=three_otu.names)
    app = get_app("model", c3_model, tree=tree_topology)
    expected = app(three_otu)

    model = Model(iq_model)

    got = piqtree.fit_tree(three_otu, tree_topology, model)
    check_likelihood(got, expected)
    check_motif_probs(got, expected.tree)
    check_rate_parameters(got, expected.tree)
    check_branch_lengths(got, expected.tree)
    check_model_name(got, str(model))


@pytest.mark.parametrize(
    ("iq_model", "c3_model"),
    [
        (StandardDnaModel.JC, "JC69"),
        (StandardDnaModel.K80, "K80"),
        (StandardDnaModel.GTR, "GTR"),
        (StandardDnaModel.TN, "TN93"),
        (StandardDnaModel.HKY, "HKY85"),
        (StandardDnaModel.F81, "F81"),
    ],
)
def test_fit_tree_str_model(
    three_otu: Alignment,
    iq_model: StandardDnaModel,
    c3_model: str,
) -> None:
    tree_topology = make_tree(tip_names=three_otu.names)
    app = get_app("model", c3_model, tree=tree_topology)
    expected = app(three_otu)

    model = str(Model(iq_model))

    got = piqtree.fit_tree(three_otu, tree_topology, model)
    check_likelihood(got, expected)
    check_motif_probs(got, expected.tree)
    check_rate_parameters(got, expected.tree)
    check_branch_lengths(got, expected.tree)


def test_fit_tree_fixed_branch_length(
    three_otu: Alignment,
) -> None:
    tree_topology: PhyloNode = make_tree(tip_names=three_otu.names)
    lengths = (0.1, 0.2, 0.3)
    for i, node in enumerate(tree_topology.postorder(include_self=False)):
        node.length = lengths[i]

    with_fixed = piqtree.fit_tree(
        three_otu,
        tree_topology,
        "GTR",
        bl_fixed=True,
    )
    without_fixed = piqtree.fit_tree(
        three_otu,
        tree_topology,
        "GTR",
        bl_fixed=False,
    )

    assert "lnL" in with_fixed.params
    assert "lnL" in without_fixed.params

    assert with_fixed.params["lnL"] != pytest.approx(without_fixed.params["lnL"])

    for fixed_node, not_fixed_node, original_node in zip(
        with_fixed.postorder(include_self=False),
        without_fixed.postorder(include_self=False),
        tree_topology.postorder(include_self=False),
        strict=True,
    ):
        assert fixed_node.name == not_fixed_node.name
        assert not_fixed_node.name == original_node.name

        assert fixed_node.length == original_node.length
        assert not_fixed_node.length != original_node.length


@pytest.mark.parametrize(
    "model_str",
    [
        "GTR{4.39,5.30,4.39,1.0,12.1}+F{0.1,0.2,0.3,0.4}+I{0.2}+G3{0.7}",
        "GTR{4.39,5.30,4.39,1.0,12.1,3.2}+R3{0.1,0.8,0.2,0.5,0.7,0.6}",
        "WS3.3b{0.5,-0.2}+F{0.6,0.1,0.2,0.1}+I{0.1}",
    ],
)
def test_fit_tree_paramaterisation(three_otu: Alignment, model_str: str) -> None:
    tree_topology = make_tree(tip_names=three_otu.names)

    tree = piqtree.fit_tree(three_otu, tree_topology, model_str)

    assert isinstance(tree.params["lnL"], float)
    for node in tree.preorder(include_self=False):
        assert cast("float", node.length) > 0


def test_special_characters(three_otu: Alignment) -> None:
    def _renamer(before: str) -> str:
        if before == three_otu.names[0]:
            return r"_F''<.'_l_?|\y}_F_o_!@#$%^&*x_''"
        return before

    three_otu = three_otu.renamed_seqs(_renamer)

    to_fit = make_tree(r"(_F''<.'_l_?|\y}_F_o_!@#$%^&*x_'', (Rhesus, Mouse))")
    tree = piqtree.fit_tree(three_otu, to_fit, "GTR")

    assert isinstance(tree.params["lnL"], float)
    for node in tree.preorder(include_self=False):
        assert cast("float", node.length) > 0

    assert set(three_otu.names) == set(tree.get_tip_names())


def test_fit_tree_other_options(four_otu: Alignment) -> None:
    tree = make_tree(tip_names=four_otu.names)
    fitted_tree = piqtree.fit_tree(
        four_otu,
        tree,
        "GTR",
        other_options="-blmin 0.1 -blmax 0.5",
    )
    for node in fitted_tree.postorder(include_self=False):
        assert node.length is not None
        assert 0.1 <= node.length <= 0.5

    with pytest.raises(IqTreeError):
        # Invalid extra options
        fitted_tree = piqtree.fit_tree(
            four_otu,
            tree,
            "GTR",
            other_options="-blminn 0.1 -blmaax 0.5",
        )
