import re
from collections.abc import Callable

import pytest

from piqtree import TreeGenMode, random_tree
from piqtree.exceptions import IqTreeError


@pytest.mark.parametrize("num_taxa", [10, 50, 100])
@pytest.mark.parametrize("tree_mode", list(TreeGenMode))
def test_random_tree(num_taxa: int, tree_mode: TreeGenMode) -> None:
    tree = random_tree(num_taxa, tree_mode, rand_seed=1)
    assert len(tree.tips()) == num_taxa


@pytest.mark.parametrize("rand_seed", [0, 1234])
def test_random_tree_determinism(rand_seed: int) -> None:
    tree_1 = random_tree(10, TreeGenMode.YULE_HARDING, rand_seed=rand_seed)
    tree_2 = random_tree(10, TreeGenMode.YULE_HARDING, rand_seed=rand_seed)

    for node_a, node_b in zip(tree_1.postorder(), tree_2.postorder(), strict=True):
        assert node_a.name == node_b.name
        assert node_a.length == node_b.length
        assert len(node_a.children) == len(node_b.children)


@pytest.mark.parametrize("num_taxa", [10, 50, 100])
@pytest.mark.parametrize("tree_mode", list(TreeGenMode))
def test_random_tree_no_seed(num_taxa: int, tree_mode: TreeGenMode) -> None:
    tree = random_tree(num_taxa, tree_mode)
    assert len(tree.tips()) == num_taxa


@pytest.mark.parametrize("num_taxa", [-1, 0, 1, 2])
@pytest.mark.parametrize("tree_mode", list(TreeGenMode))
def test_invalid_taxa(
    num_taxa: int,
    tree_mode: TreeGenMode,
) -> None:
    with pytest.raises(IqTreeError):
        _ = random_tree(num_taxa, tree_mode, rand_seed=1)


@pytest.mark.parametrize("tree_mode", list(TreeGenMode))
@pytest.mark.parametrize("transform", [str.upper, str.lower])
def test_random_tree_str_mode(
    tree_mode: TreeGenMode,
    transform: Callable[[str], str],
) -> None:
    expected = random_tree(10, tree_mode, rand_seed=1)
    got = random_tree(10, transform(tree_mode.name), rand_seed=1)
    assert str(got) == str(expected)


@pytest.mark.parametrize("tree_mode", ["YULE", "", "BIRTH-DEATH"])
def test_random_tree_bad_str_mode(tree_mode: str) -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(f"Unknown tree generation mode {tree_mode!r}"),
    ):
        _ = random_tree(10, tree_mode, rand_seed=1)
