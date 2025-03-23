"""Test combinations of calls which under previous versions resulted in a segmentation fault."""

import pytest
from cogent3 import make_aligned_seqs, make_tree

from piqtree import TreeGenMode, build_tree, fit_tree, random_trees
from piqtree.exceptions import IqTreeError
from piqtree.model import DiscreteGammaModel, DnaModel, FreeRateModel, Model


def test_two_build_random_trees() -> None:
    """
    Calling build tree twice followed by random trees with a bad input
    used to result in a Segmentation Fault in a previous version.
    """
    aln = make_aligned_seqs({"a": "GGG", "b": "GGC", "c": "AAC", "d": "AAA"})

    build_tree(aln, Model(DnaModel.JC), 1)
    build_tree(aln, Model(DnaModel.JC), 2)

    with pytest.raises(IqTreeError):
        random_trees(3, 2, TreeGenMode.BALANCED, 1)


def test_two_fit_random_trees() -> None:
    """
    Calling fit tree twice followed by random trees with a bad input
    used to result in a Segmentation Fault in a previous version.
    """
    aln = make_aligned_seqs({"a": "GGG", "b": "GGC", "c": "AAC", "d": "AAA"})
    tree = make_tree("(a,b,(c,d));")

    fit_tree(aln, tree, Model(DnaModel.JC), 1)
    fit_tree(aln, tree, Model(DnaModel.JC), 2)

    with pytest.raises(IqTreeError):
        random_trees(3, 2, TreeGenMode.BALANCED, 1)


@pytest.mark.parametrize("rate_model_class", [DiscreteGammaModel, FreeRateModel])
@pytest.mark.parametrize("categories", [0, -4])
def test_two_invalid_models(
    rate_model_class: type[DiscreteGammaModel] | type[FreeRateModel],
    categories: int,
) -> None:
    """
    Calling build_tree multiple times with an invalid
    model has resulted in a Segmentation Fault.
    """
    aln = make_aligned_seqs({"a": "GGG", "b": "GGC", "c": "AAC", "d": "AAA"})

    with pytest.raises(IqTreeError):
        _ = build_tree(aln, Model(DnaModel.JC, rate_model=rate_model_class(categories)))

    with pytest.raises(IqTreeError):
        _ = build_tree(aln, Model(DnaModel.JC, rate_model=rate_model_class(categories)))
