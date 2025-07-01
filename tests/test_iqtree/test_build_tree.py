import re

import pytest
from cogent3 import make_tree
from cogent3.core.alignment import Alignment

import piqtree
from piqtree.exceptions import IqTreeError
from piqtree.model import (
    DiscreteGammaModel,
    FreeRateModel,
    FreqType,
    LieModel,
    LieModelInstance,
    Model,
    RateModel,
    StandardDnaModel,
    SubstitutionModel,
)


def check_build_tree(
    aln: Alignment,
    dna_model: SubstitutionModel,
    freq_type: FreqType | None = None,
    rate_model: RateModel | None = None,
    *,
    invariable_sites: bool = False,
    coerce_str: bool = False,
) -> None:
    expected = make_tree("(Human,Chimpanzee,(SpermWhale,HumpbackW));")

    model = Model(
        dna_model,
        freq_type=freq_type if freq_type else None,
        invariable_sites=invariable_sites,
        rate_model=rate_model,
    )

    got = piqtree.build_tree(aln, str(model) if coerce_str else model)
    # Check topology
    assert expected.same_topology(got.unrooted())
    # Check if branch lengths exist
    assert all("length" in v.params for v in got.get_edge_vector())


@pytest.mark.parametrize("dna_model", StandardDnaModel.iter_available_models())
@pytest.mark.parametrize("freq_type", list(FreqType))
def test_non_lie_build_tree(
    four_otu: Alignment,
    dna_model: StandardDnaModel,
    freq_type: FreqType,
) -> None:
    check_build_tree(four_otu, dna_model, freq_type)


@pytest.mark.parametrize("lie_model", LieModel.iter_available_models())
def test_lie_build_tree(four_otu: Alignment, lie_model: LieModelInstance) -> None:
    check_build_tree(four_otu, lie_model)


@pytest.mark.parametrize("lie_model", LieModel.iter_available_models()[:3])
def test_str_build_tree(four_otu: Alignment, lie_model: LieModelInstance) -> None:
    check_build_tree(four_otu, lie_model, coerce_str=True)


@pytest.mark.parametrize("dna_model", StandardDnaModel.iter_available_models()[:3])
@pytest.mark.parametrize("invariable_sites", [False, True])
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
    four_otu: Alignment,
    dna_model: StandardDnaModel,
    invariable_sites: bool,
    rate_model: RateModel,
) -> None:
    check_build_tree(
        four_otu,
        dna_model,
        rate_model=rate_model,
        invariable_sites=invariable_sites,
    )


def test_build_tree_inadequate_bootstrapping(four_otu: Alignment) -> None:
    with pytest.raises(IqTreeError, match=re.escape("#replicates must be >= 1000")):
        piqtree.build_tree(
            four_otu,
            Model(StandardDnaModel.GTR),
            bootstrap_replicates=10,
        )


def test_build_tree_bootstrapping(four_otu: Alignment) -> None:
    tree = piqtree.build_tree(
        four_otu,
        Model(StandardDnaModel.GTR),
        bootstrap_replicates=1000,
    )

    supported_node = max(tree.children, key=lambda x: len(x.children))
    assert "support" in supported_node.params
