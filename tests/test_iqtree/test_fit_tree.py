import numpy as np
import pytest
from cogent3 import get_app, make_tree
from cogent3.app.result import model_result
from cogent3.core.alignment import Alignment
from cogent3.core.tree import PhyloNode

import piqtree
from piqtree.model import Model, StandardDnaModel


def check_likelihood(got: PhyloNode, expected: model_result) -> None:
    assert got.params["lnL"] == pytest.approx(expected.lnL)


def check_motif_probs(got: PhyloNode, expected: PhyloNode) -> None:
    expected = expected.params["mprobs"]
    got = got.params["mprobs"]

    expected_keys = set(expected.keys())
    got_keys = set(got.keys())

    # Check that the base characters are the same
    assert expected_keys == got_keys

    # Check that the probs are the same
    expected_values = [expected[key] for key in expected_keys]
    got_values = [got[key] for key in expected_keys]
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
    got = got.tip_to_tip_distances()
    expected = expected.tip_to_tip_distances()
    # make sure the distance matrices have the same name order
    # so we can just compare entire numpy arrays
    expected = expected.take_dists(got.names)
    # Check that the keys of branch lengths are the same
    assert set(got.names) == set(expected.names)

    # Check that the branch lengths are the same
    np.testing.assert_allclose(got.array, expected.array, atol=1e-4)


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
        assert node.length > 0


def test_special_characters(three_otu: Alignment) -> None:
    def _renamer(before: str) -> str:
        if before == three_otu.names[0]:
            return r"_F''<.'_l_?|\y}_F_o_!@#$%^&*x_''"
        return before

    three_otu = three_otu.rename_seqs(_renamer)

    to_fit = make_tree(r"(_F''<.'_l_?|\y}_F_o_!@#$%^&*x_'', (Rhesus, Mouse))")
    tree = piqtree.fit_tree(three_otu, to_fit, "GTR")

    assert isinstance(tree.params["lnL"], float)
    for node in tree.preorder(include_self=False):
        assert node.length > 0

    assert set(three_otu.names) == set(tree.get_tip_names())
