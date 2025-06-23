import pytest
from cogent3.core.alignment import Alignment

import piqtree
from piqtree.model import (
    AaModel,
    Model,
    SubstitutionModel,
)


def check_build_tree(
    aln: Alignment,
    aa_model: SubstitutionModel,
    *,
    coerce_str: bool = False,
) -> None:
    model = Model(aa_model)

    got = piqtree.build_tree(aln, str(model) if coerce_str else model, rand_seed=1)

    # Check if all branch lengths exist
    assert all("length" in v.params for v in got.get_edge_vector())


@pytest.mark.parametrize("aa_model", AaModel.iter_available_models())
def test_protein_build_tree(protein_four_otu: Alignment, aa_model: AaModel) -> None:
    check_build_tree(protein_four_otu, aa_model, coerce_str=True)


@pytest.mark.parametrize("aa_model", AaModel.iter_available_models())
def test_str_protein_build_tree(protein_four_otu: Alignment, aa_model: AaModel) -> None:
    check_build_tree(protein_four_otu, aa_model, coerce_str=False)
