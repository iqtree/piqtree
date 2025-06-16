import re

import numpy as np
import pytest
from cogent3 import make_tree
from cogent3.core.alignment import Alignment

from piqtree import jc_distances, nj_tree


def test_nj_tree(five_otu: Alignment) -> None:
    expected = make_tree("(((Human, Chimpanzee), Rhesus), Manatee, Dugong);")

    dists = jc_distances(five_otu)
    actual = nj_tree(dists)

    assert expected.same_topology(actual)


def test_nj_tree_allow_negative(all_otu: Alignment) -> None:
    # a distance matrix can produce trees with negative branch lengths
    dists = jc_distances(all_otu)

    # check that all branch lengths are non-negative, by default
    tree1 = nj_tree(dists)
    assert all(node.length >= 0 for node in tree1.traverse(include_self=False))

    # check that some branch lengths are negative when allow_negative=True
    tree2 = nj_tree(dists, allow_negative=True)
    assert any(node.length < 0 for node in tree2.traverse(include_self=False))


def test_nj_tree_nan(four_otu: Alignment) -> None:
    dists = jc_distances(four_otu)
    dists[1, 0] = dists[0, 1] = np.nan

    with pytest.raises(
        ValueError,
        match=re.escape("The pairwise distance matrix cannot contain NaN values."),
    ):
        nj_tree(dists, allow_negative=True)
