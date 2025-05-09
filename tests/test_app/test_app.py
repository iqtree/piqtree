import pytest
from cogent3 import get_app, make_tree
from cogent3.core.new_alignment import Alignment

import piqtree
from piqtree import jc_distances, make_model


def test_piqtree_phylo(four_otu: Alignment) -> None:
    expected = make_tree("(Human,Chimpanzee,(Rhesus,Mouse));")
    app = get_app("piqtree_phylo", model="JC")
    got = app(four_otu)
    assert expected.same_topology(got)


def test_piqtree_phylo_support(four_otu: Alignment) -> None:
    app = get_app("piqtree_phylo", model=make_model("JC"), bootstrap_reps=1000)
    got = app(four_otu)
    supports = [
        node.params.get("support", None)
        for node in got.postorder()
        if not node.is_tip() and node.name != "root"
    ]
    assert all(supports)


def test_piqtree_fit(three_otu: Alignment) -> None:
    tree = make_tree(tip_names=three_otu.names)
    app = get_app("model", "JC69", tree=tree)
    expected = app(three_otu)
    piphylo = get_app("piqtree_fit", tree=tree, model="JC")
    got = piphylo(three_otu)
    assert got.params["lnL"] == pytest.approx(expected.lnL)


@pytest.mark.parametrize("num_trees", [1, 10, 20])
@pytest.mark.parametrize("num_taxa", [10, 50, 100])
@pytest.mark.parametrize("tree_mode", list(piqtree.TreeGenMode))
def test_piqtree_random_trees(
    num_trees: int,
    num_taxa: int,
    tree_mode: piqtree.TreeGenMode,
) -> None:
    app = get_app(
        "piqtree_random_trees",
        tree_mode=tree_mode,
        num_trees=num_trees,
        rand_seed=1,
    )
    trees = app(num_taxa)
    assert len(trees) == num_trees

    for tree in trees:
        assert len(tree.tips()) == num_taxa


def test_piqtree_jc_distances(five_otu: Alignment) -> None:
    app = get_app("piqtree_jc_dists")
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


def test_piqtree_nj(five_otu: Alignment) -> None:
    dists = jc_distances(five_otu)

    expected = make_tree("(((Human, Chimpanzee), Rhesus), Manatee, Dugong);")

    app = get_app("piqtree_nj")

    actual = app(dists)

    assert expected.same_topology(actual)


def test_mfinder(five_otu: Alignment) -> None:
    from piqtree.iqtree import ModelFinderResult

    app = get_app("piqtree_mfinder")
    got = app(five_otu)
    assert isinstance(got, ModelFinderResult)


def test_mfinder_result_roundtrip(five_otu: Alignment) -> None:
    from piqtree.iqtree import ModelFinderResult

    app = get_app("piqtree_mfinder")
    got = app(five_otu)
    rd = got.to_rich_dict()
    inflated = ModelFinderResult.from_rich_dict(rd)
    assert isinstance(inflated, ModelFinderResult)
    assert str(got.best_aicc) == str(inflated.best_aicc)


@pytest.mark.parametrize("use_hook", [None, "piqtree"])
def test_quick_tree_hook(four_otu: Alignment, use_hook: str | None) -> None:
    tree = four_otu.quick_tree(use_hook=use_hook)
    assert tree.params["provenance"] == "piqtree"
