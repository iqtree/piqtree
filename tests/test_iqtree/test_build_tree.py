import re

import pytest
from cogent3 import ArrayAlignment, make_tree

import piqtree
from piqtree.exceptions import IqTreeError
from piqtree.model import (
    DiscreteGammaModel,
    DnaModel,
    FreeRateModel,
    FreqType,
    Model,
    RateModel,
)


def check_build_tree(
    four_otu: ArrayAlignment,
    dna_model: DnaModel,
    freq_type: FreqType | None = None,
    rate_model: RateModel | None = None,
    *,
    invariant_sites: bool = False,
) -> None:
    expected = make_tree("(Human,Chimpanzee,(Rhesus,Mouse));")

    model = Model(
        dna_model,
        freq_type=freq_type if freq_type else None,
        invariant_sites=invariant_sites,
        rate_model=rate_model,
    )

    got1 = piqtree.build_tree(four_otu, model, rand_seed=1)
    got1 = got1.unrooted()
    # Check topology
    assert expected.same_topology(got1.unrooted())
    # Check if branch lengths exist
    assert all("length" in v.params for v in got1.get_edge_vector())

    # Should be similar for any seed
    got2 = piqtree.build_tree(four_otu, model, rand_seed=None)
    got2 = got2.unrooted()
    assert expected.same_topology(got2)
    assert all("length" in v.params for v in got2.get_edge_vector())


@pytest.mark.parametrize("dna_model", list(DnaModel)[:22])
@pytest.mark.parametrize("freq_type", list(FreqType))
def test_non_lie_build_tree(
    four_otu: ArrayAlignment,
    dna_model: DnaModel,
    freq_type: FreqType,
) -> None:
    check_build_tree(four_otu, dna_model, freq_type)


@pytest.mark.parametrize("dna_model", list(DnaModel)[22:])
def test_lie_build_tree(four_otu: ArrayAlignment, dna_model: DnaModel) -> None:
    check_build_tree(four_otu, dna_model)


@pytest.mark.parametrize("dna_model", list(DnaModel)[:5])
@pytest.mark.parametrize("invariant_sites", [False, True])
@pytest.mark.parametrize(
    "rate_model",
    [
        None,
        DiscreteGammaModel(),
        FreeRateModel(),
        DiscreteGammaModel(6),
        FreeRateModel(6),
    ],
)
def test_rate_model_build_tree(
    four_otu: ArrayAlignment,
    dna_model: DnaModel,
    invariant_sites: bool,
    rate_model: RateModel,
) -> None:
    check_build_tree(
        four_otu,
        dna_model,
        rate_model=rate_model,
        invariant_sites=invariant_sites,
    )


def test_build_tree_inadequate_bootstrapping(four_otu: ArrayAlignment) -> None:
    with pytest.raises(IqTreeError, match=re.escape("#replicates must be >= 1000")):
        piqtree.build_tree(four_otu, Model(DnaModel.GTR), bootstrap_replicates=10)


def test_build_tree_bootstrapping(four_otu: ArrayAlignment) -> None:
    tree = piqtree.build_tree(four_otu, Model(DnaModel.GTR), bootstrap_replicates=1000)

    supported_node = max(tree.children, key=lambda x: len(x.children))
    assert "support" in supported_node.params
