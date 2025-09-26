import numpy as np
from cogent3 import make_tree
from numpy.testing import assert_array_equal

import piqtree


def test_robinson_foulds() -> None:
    tree1 = make_tree("(A,B,(C,D));")
    tree2 = make_tree("(A,C,(B,D));")
    pairwise_distances = piqtree.robinson_foulds([tree1, tree2])
    assert_array_equal(pairwise_distances, np.array([[0, 2], [2, 0]]))
