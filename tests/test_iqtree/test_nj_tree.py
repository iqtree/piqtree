import pathlib

from cogent3 import ArrayAlignment, make_tree
from cogent3.util.deserialise import deserialise_object

from piqtree import jc_distances, nj_tree


def test_nj_tree(five_otu: ArrayAlignment) -> None:
    expected = make_tree("(((Human, Chimpanzee), Rhesus), Manatee, Dugong);")

    dists = jc_distances(five_otu)
    actual = nj_tree(dists)

    assert expected.same_topology(actual)


def test_nj_tree_allow_negative(DATA_DIR: pathlib.Path) -> None:
    # a distance matrix can produce trees with negative branch lengths
    dists = deserialise_object(DATA_DIR / "distance_matrix.json")

    # check that all branch lengths are non-negative, by default
    tree1 = nj_tree(dists)
    lengths1 = [v.params["length"] for v in tree1.get_edge_vector() if v.name != "root"]
    assert all(length >= 0 for length in lengths1)

    # check that some branch lengths are negative when allow_negative=True
    tree2 = nj_tree(dists, allow_negative=True)
    lengths2 = [v.params["length"] for v in tree2.get_edge_vector() if v.name != "root"]
    assert any(length < 0 for length in lengths2)
