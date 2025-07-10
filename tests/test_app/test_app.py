import pytest
from cogent3 import PhyloNode, get_app, make_tree
from cogent3.core.alignment import Alignment

import piqtree
from piqtree import ModelFinderResult, jc_distances, make_model


def test_piq_build_tree(four_otu: Alignment) -> None:
    expected = make_tree("(Human,Chimpanzee,(SpermWhale,HumpbackW));")
    app = get_app("piq_build_tree", model="JC")
    got = app(four_otu)
    assert expected.same_topology(got)
    assert got.source == four_otu.source


def test_piq_build_tree_support(four_otu: Alignment) -> None:
    app = get_app("piq_build_tree", model=make_model("JC"), bootstrap_reps=1000)
    got = app(four_otu)
    supports = [
        node.params.get("support", None)
        for node in got.postorder()
        if not node.is_tip() and node.name != "root"
    ]
    assert all(supports)


def test_piq_fit_tree(three_otu: Alignment) -> None:
    tree = make_tree(tip_names=three_otu.names)
    app = get_app("model", "JC69", tree=tree)
    expected = app(three_otu)
    piphylo = get_app("piq_fit_tree", tree=tree, model="JC")
    got = piphylo(three_otu)
    assert got.params["lnL"] == pytest.approx(expected.lnL)
    assert got.source == three_otu.source


@pytest.mark.parametrize("num_taxa", [10, 50, 100])
@pytest.mark.parametrize("tree_mode", list(piqtree.TreeGenMode))
def test_piq_random_tree(num_taxa: int, tree_mode: piqtree.TreeGenMode) -> None:
    app = get_app("piq_random_tree", tree_mode=tree_mode, rand_seed=1)

    tree = app(num_taxa)
    assert len(tree.tips()) == num_taxa


def test_piq_jc_distances(five_otu: Alignment) -> None:
    app = get_app("piq_jc_distances")
    dists = app(five_otu)

    assert (
        0 < dists["Human", "Chimpanzee"] < dists["Human", "Dugong"]
    )  # chimpanzee closer than rhesus
    assert (
        0 < dists["Human", "Rhesus"] < dists["Human", "Manatee"]
    )  # rhesus closer than manatee
    assert (
        0 < dists["Human", "Rhesus"] < dists["Human", "Dugong"]
    )  # rhesus closer than dugong

    assert (
        0 < dists["Chimpanzee", "Rhesus"] < dists["Chimpanzee", "Manatee"]
    )  # rhesus closer than manatee
    assert (
        0 < dists["Chimpanzee", "Rhesus"] < dists["Chimpanzee", "Dugong"]
    )  # rhesus closer than dugong

    assert (
        0 < dists["Manatee", "Dugong"] < dists["Manatee", "Rhesus"]
    )  # dugong closer than rhesus

    assert dists.source == five_otu.source


def test_piq_nj_tree(five_otu: Alignment) -> None:
    dists = jc_distances(five_otu)

    expected = make_tree("(((Human, Chimpanzee), Rhesus), Manatee, Dugong);")

    app = get_app("piq_nj_tree")

    actual = app(dists)

    assert expected.same_topology(actual)


def test_piq_model_finder(five_otu: Alignment) -> None:
    app = get_app("piq_model_finder")
    got = app(five_otu)
    assert isinstance(got, ModelFinderResult)


def test_piq_model_finder_result_roundtrip(five_otu: Alignment) -> None:
    app = get_app("piq_model_finder")
    got = app(five_otu)
    rd = got.to_rich_dict()
    inflated = ModelFinderResult.from_rich_dict(rd)
    assert isinstance(inflated, ModelFinderResult)
    assert str(got.best_aicc) == str(inflated.best_aicc)


def tree_equal(tree1: PhyloNode, tree2: PhyloNode) -> bool:
    return str(tree1.sorted()) == str(tree2.sorted())


def test_piq_consesus_tree(five_trees: list[PhyloNode]) -> None:
    app_majority = get_app("piq_consensus_tree")
    app_strict = get_app("piq_consensus_tree", min_support=1)
    app_0_3 = get_app("piq_consensus_tree", min_support=0.3)

    expected_majority = make_tree("((a,b),(((e,f),d),c))")
    got_majority = app_majority(five_trees)
    assert tree_equal(got_majority, expected_majority)

    expected_strict = make_tree("(a,b,c,d,(e,f));")
    got_strict = app_strict(five_trees)
    assert tree_equal(got_strict, expected_strict)

    expected_0_3 = make_tree("((a,b),(((e,f),d),c))")
    got_0_3 = app_0_3(five_trees)
    assert tree_equal(got_0_3, expected_0_3)


def test_quick_tree_hook(four_otu: Alignment) -> None:
    tree = four_otu.quick_tree(use_hook="piqtree")
    assert tree.params["provenance"] == "piqtree"
