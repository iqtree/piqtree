import re

import pytest
from cogent3 import make_tree
from cogent3.core.alignment import Alignment

import piqtree
from piqtree import make_model
from piqtree.exceptions import IqTreeError
from piqtree.model import (
    CustomBaseFreq,
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


def check_build_tree_model(
    aln: Alignment,
    model: Model,
    *,
    coerce_str: bool = False,
) -> None:
    expected = make_tree("(Human,Chimpanzee,(SpermWhale,HumpbackW));")

    got = piqtree.build_tree(aln, str(model) if coerce_str else model)
    # Check topology
    assert expected.same_topology(got.unrooted())
    # Check if branch lengths exist
    assert all(v.length is not None for v in got.get_edge_vector(include_root=False))
    assert got.params["model"] == str(model)


def check_build_tree(
    aln: Alignment,
    dna_model: SubstitutionModel,
    freq_type: FreqType | CustomBaseFreq | None = None,
    rate_model: RateModel | None = None,
    *,
    invariable_sites: bool = False,
    coerce_str: bool = False,
) -> None:
    model = Model(
        dna_model,
        freq_type=freq_type if freq_type else None,
        invariable_sites=invariable_sites,
        rate_model=rate_model,
    )

    check_build_tree_model(aln, model, coerce_str=coerce_str)


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
    assert supported_node.support is not None


@pytest.mark.parametrize(
    "model_str",
    [
        "GTR{4.39,5.30,4.39,1.0,12.1}+F{0.1,0.2,0.3,0.4}+I{0.2}+G3{0.7}",
        "GTR{4.39,5.30,4.39,1.0,12.1,3.2}+R3{0.1,0.8,0.2,0.5,0.7,0.6}",
        "WS3.3b{0.5,-0.2}+F{0.6,0.1,0.2,0.1}+I{0.1}",
    ],
)
def test_build_tree_paramaterisation(four_otu: Alignment, model_str: str) -> None:
    model = make_model(model_str)
    check_build_tree_model(four_otu, model)


def test_invalid_protein_base_freq(four_otu: Alignment) -> None:
    base_freqs = ["0.05"] * 20

    model = f"K81+F{{{','.join(base_freqs)}}}"

    with pytest.raises(IqTreeError):
        _ = piqtree.build_tree(four_otu, model)


def test_too_many_dna_params(four_otu: Alignment) -> None:
    with pytest.raises(IqTreeError):
        _ = piqtree.build_tree(four_otu, "GTR{4.39,5.30,4.39,1.0,12.1,3.2,1.5}")


def test_too_few_dna_params(four_otu: Alignment) -> None:
    with pytest.raises(IqTreeError):
        _ = piqtree.build_tree(four_otu, "GTR{4.39,5.30,4.39,1.0}")


def test_build_tree_other_options(four_otu: Alignment) -> None:
    _ = piqtree.build_tree(four_otu, "GTR", other_options="--epsilon 0.1 -n 0")

    with pytest.raises(IqTreeError):
        # Invalid extra options
        _ = piqtree.build_tree(four_otu, "GTR", other_options="--epslon 0.1 -n 0")
