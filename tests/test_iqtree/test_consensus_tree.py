import re
from collections.abc import Iterable

import pytest
from cogent3 import PhyloNode, make_tree

from piqtree import consensus_tree
from piqtree.exceptions import IqTreeError


def tree_equal(tree1: PhyloNode, tree2: PhyloNode) -> bool:
    return str(tree1.sorted()) == str(tree2.sorted())


@pytest.fixture
def standard_trees() -> list[PhyloNode]:
    tree1 = make_tree("(a,(b,(c,(d,(e,f)))))")
    tree2 = make_tree("(a,(b,(c,(d,(e,f)))))")
    tree3 = make_tree("((a,b),(c,(d,(e,f))))")
    tree4 = make_tree("(((a,b),c),(d,(e,f)))")
    tree5 = make_tree("((((a,b),c),d),(e,f))")

    return [tree1, tree2, tree3, tree4, tree5]


def test_majority_consensus_tree(standard_trees: list[PhyloNode]) -> None:
    expected = make_tree("((a,b),(((e,f),d),c))")
    got_default = consensus_tree(standard_trees)
    assert tree_equal(got_default, expected)

    got = consensus_tree(standard_trees, min_support=0.5)
    assert tree_equal(got, expected)


def test_higher_support(standard_trees: list[PhyloNode]) -> None:
    expected_0_7 = make_tree("(a,b,c,(d,(e,f)));")
    got_0_7 = consensus_tree(standard_trees, min_support=0.7)
    assert tree_equal(got_0_7, expected_0_7)

    expected_0_9 = make_tree("(a,b,c,d,(e,f));")
    got_0_9 = consensus_tree(standard_trees, min_support=0.9)
    assert tree_equal(got_0_9, expected_0_9)


def test_strict_consensus_tree(standard_trees: list[PhyloNode]) -> None:
    expected = make_tree("(a,b,c,d,(e,f));")
    got = consensus_tree(standard_trees, min_support=1)
    assert tree_equal(got, expected)


def test_extended_majority_rule() -> None:
    tree1 = make_tree("(a,(b,(c,(d,e))))")
    tree2 = make_tree("(a,((b,c),(d,e)))")
    tree3 = make_tree("(a,((b,c),(d,e)))")
    tree4 = make_tree("(a,(c,(b,(d,e))))")
    tree5 = make_tree("(c,(b,(a,(d,e))))")

    got = consensus_tree([tree1, tree2, tree3, tree4, tree5], min_support=0)
    expected = make_tree("(a,((b,c),(d,e)))")

    assert tree_equal(got, expected)


def test_lower_support(standard_trees: list[PhyloNode]) -> None:
    expected = make_tree("((a,b),(((e,f),d),c))")
    for support in 0.1, 0.3:
        got = consensus_tree(standard_trees, min_support=support)
        assert tree_equal(got, expected)


def test_single_tree(standard_trees: list[PhyloNode]) -> None:
    single_tree = standard_trees[0]
    got = consensus_tree([single_tree])
    assert tree_equal(got, single_tree)


def test_even_majority_rule() -> None:
    tree1 = make_tree("(a,(b,(c,d)))")
    tree2 = make_tree("(a,(b,(c,d)))")
    tree3 = make_tree("(b,(a,(c,d)))")
    tree4 = make_tree("(c,(a,(b,d)))")

    expected = make_tree("(a,b,(c,d))")  # Not including (b, (c,d)) as that's only in 2

    got = consensus_tree([tree1, tree2, tree3, tree4])

    assert tree_equal(got, expected)


@pytest.mark.parametrize("min_support", [-1, -0.1, 1.00001, 2.5])
def test_bad_min_support(standard_trees: list[PhyloNode], min_support: float) -> None:
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Only min support values in the range 0 <= value < 1 are supported, got {min_support}",
        ),
    ):
        consensus_tree(standard_trees, min_support=min_support)


@pytest.mark.parametrize(
    "trees",
    [
        [make_tree("(a,(b,(c,d)))"), make_tree("(a,(b,(c,(d,e))))")],
        [make_tree("(a,(b,(c,(d,e))))"), make_tree("(a,(b,(c,d)))")],
        (make_tree("(a,(b,(c,e)))"), make_tree("(a,(b,(c,d)))")),
        (make_tree("(a,(b,(c,d)))"), make_tree("(a,(b,(c,a)))")),
    ],
)
def test_bad_trees(trees: Iterable[PhyloNode]) -> None:
    with pytest.raises(ValueError, match=re.escape("Trees must be on same taxa set.")):
        consensus_tree(trees)


def test_no_trees() -> None:
    with pytest.raises(IqTreeError):
        consensus_tree([])
