import pytest
from cogent3.core.alignment import Alignment

import piqtree
from piqtree import make_model
from piqtree.exceptions import IqTreeError
from piqtree.model import (
    AaModel,
    Model,
    SubstitutionModel,
)


def check_build_tree_model(
    aln: Alignment,
    model: Model,
    *,
    coerce_str: bool = False,
) -> None:
    got = piqtree.build_tree(aln, str(model) if coerce_str else model, rand_seed=1)

    # Check if all branch lengths exist
    assert all("length" in v.params for v in got.get_edge_vector())


def check_build_tree(
    aln: Alignment,
    aa_model: SubstitutionModel,
    *,
    coerce_str: bool = False,
) -> None:
    model = Model(aa_model)
    check_build_tree_model(aln, model, coerce_str=coerce_str)


@pytest.mark.parametrize("aa_model", AaModel.iter_available_models())
def test_protein_build_tree(protein_four_otu: Alignment, aa_model: AaModel) -> None:
    check_build_tree(protein_four_otu, aa_model, coerce_str=True)


@pytest.mark.parametrize("aa_model", AaModel.iter_available_models())
def test_str_protein_build_tree(protein_four_otu: Alignment, aa_model: AaModel) -> None:
    check_build_tree(protein_four_otu, aa_model, coerce_str=False)


def test_parameterised_model(protein_four_otu: Alignment) -> None:
    base_freq_values = ["0.05"] * 20

    model = make_model(f"rtREV+F{{{','.join(base_freq_values)}}}+I{{0.2}}+G3{{1.2}}")
    check_build_tree_model(protein_four_otu, model)


def test_invalid_protein_base_freq(protein_four_otu: Alignment) -> None:
    model = "FLAVI+F{0.2,0.3,0.4,0.1}"

    with pytest.raises(IqTreeError):
        _ = piqtree.build_tree(protein_four_otu, model)
