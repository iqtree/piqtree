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
    got = got.get_distances()
    expected = expected.get_distances()
    # Check that the keys of branch lengths are the same
    assert got.keys() == expected.keys()

    # Check that the branch lengths are the same
    expected_values = [expected[key] for key in expected]
    got_values = [got[key] for key in expected]

    assert all(
        got == pytest.approx(exp, rel=1e-2)
        for got, exp in zip(got_values, expected_values, strict=True)
    )


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

    got1 = piqtree.fit_tree(three_otu, tree_topology, model)
    check_likelihood(got1, expected)
    check_motif_probs(got1, expected.tree)
    check_rate_parameters(got1, expected.tree)
    check_branch_lengths(got1, expected.tree)


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

    got1 = piqtree.fit_tree(three_otu, tree_topology, model)
    check_likelihood(got1, expected)
    check_motif_probs(got1, expected.tree)
    check_rate_parameters(got1, expected.tree)
    check_branch_lengths(got1, expected.tree)
