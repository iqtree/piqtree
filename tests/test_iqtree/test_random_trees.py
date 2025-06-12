import pytest

import piqtree
import piqtree.exceptions


@pytest.mark.parametrize("num_taxa", [10, 50, 100])
@pytest.mark.parametrize("tree_mode", list(piqtree.TreeGenMode))
def test_random_tree(num_taxa: int, tree_mode: piqtree.TreeGenMode) -> None:
    tree = piqtree.random_tree(num_taxa, tree_mode, rand_seed=1)
    assert len(tree.tips()) == num_taxa


@pytest.mark.parametrize("num_taxa", [10, 50, 100])
@pytest.mark.parametrize("tree_mode", list(piqtree.TreeGenMode))
def test_random_tree_no_seed(num_taxa: int, tree_mode: piqtree.TreeGenMode) -> None:
    tree = piqtree.random_tree(num_taxa, tree_mode)
    assert len(tree.tips()) == num_taxa


@pytest.mark.parametrize("num_taxa", [-1, 0, 1, 2])
@pytest.mark.parametrize("tree_mode", list(piqtree.TreeGenMode))
def test_invalid_taxa(
    num_taxa: int,
    tree_mode: piqtree.TreeGenMode,
) -> None:
    with pytest.raises(piqtree.exceptions.IqTreeError):
        _ = piqtree.random_tree(num_taxa, tree_mode, rand_seed=1)
