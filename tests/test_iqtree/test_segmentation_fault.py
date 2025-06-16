"""Test combinations of calls which under previous versions resulted in a segmentation fault."""

import pytest
from cogent3 import make_aligned_seqs, make_tree
from cogent3.core.alignment import Alignment

from piqtree import TreeGenMode, build_tree, fit_tree, random_tree
from piqtree.exceptions import IqTreeError
from piqtree.model import DiscreteGammaModel, DnaModel, FreeRateModel, Model


@pytest.fixture
def tiny_alignment() -> Alignment:
    return make_aligned_seqs(
        {"a": "GGG", "b": "GGC", "c": "AAC", "d": "AAA"},
        moltype="dna",
    )


def test_two_build_random_tree(tiny_alignment: Alignment) -> None:
    """
    Calling build tree twice followed by random tree with a bad input
    used to result in a Segmentation Fault in a previous version.
    """

    build_tree(tiny_alignment, Model(DnaModel.JC), 1)
    build_tree(tiny_alignment, Model(DnaModel.JC), 2)

    with pytest.raises(IqTreeError):
        random_tree(2, TreeGenMode.BALANCED, 1)


def test_two_fit_random_tree(tiny_alignment: Alignment) -> None:
    """
    Calling fit tree twice followed by random tree with a bad input
    used to result in a Segmentation Fault in a previous version.
    """
    tree = make_tree("(a,b,(c,d));")

    fit_tree(tiny_alignment, tree, Model(DnaModel.JC), 1)
    fit_tree(tiny_alignment, tree, Model(DnaModel.JC), 2)

    with pytest.raises(IqTreeError):
        random_tree(2, TreeGenMode.BALANCED, 1)


@pytest.mark.parametrize("rate_model_class", [DiscreteGammaModel, FreeRateModel])
@pytest.mark.parametrize("categories", [0, -4])
def test_two_invalid_models(
    tiny_alignment: Alignment,
    rate_model_class: type[DiscreteGammaModel] | type[FreeRateModel],
    categories: int,
) -> None:
    """
    Calling build_tree multiple times with an invalid
    model has resulted in a Segmentation Fault.
    """
    with pytest.raises(IqTreeError):
        _ = build_tree(
            tiny_alignment,
            Model(DnaModel.JC, rate_model=rate_model_class(categories)),
        )

    with pytest.raises(IqTreeError):
        _ = build_tree(
            tiny_alignment,
            Model(DnaModel.JC, rate_model=rate_model_class(categories)),
        )
