import pytest
from cogent3 import make_tree
from cogent3.core.tree import PhyloNode

from piqtree import Model, simulate_alignment
from piqtree.model import (
    AaModel,
    LieModel,
    StandardDnaModel,
    SubstitutionModel,
)
from piqtree.model._substitution_model import LieModelInstance


@pytest.fixture
def four_taxon_unrooted_tree() -> PhyloNode:  # Unrooted
    return make_tree("((a:0.1,b:0.2):0.05,c:0.3,d:0.1);")


@pytest.fixture
def five_taxon_rooted_tree() -> PhyloNode:  # Rooted
    return make_tree("(((a:0.1,b:0.2):0.05,(c:0.3,d:0.1):0.2):0.05,e:0.4);")


def check_simulate_alignment(
    tree: PhyloNode,
    model: Model | str,
    length: int | None = None,
    insertion_rate: float = 0.0,
    deletion_rate: float = 0.0,
) -> None:
    if isinstance(model, Model):
        submod = model.submod_type
        if isinstance(submod, LieModelInstance) and submod.lie_model in (
            LieModel.LIE_1_1,
            LieModel.LIE_3_3a,
            LieModel.LIE_4_4a,
            LieModel.LIE_6_7a,
            LieModel.LIE_9_20a,
            LieModel.LIE_9_20b,
            LieModel.LIE_12_12,
        ):
            return
    if length is None:
        aln = simulate_alignment(
            tree,
            model,
            insertion_rate=insertion_rate,
            deletion_rate=deletion_rate,
        )
    else:
        aln = simulate_alignment(
            tree,
            model,
            length,
            insertion_rate=insertion_rate,
            deletion_rate=deletion_rate,
        )
    if length is None:
        length = 1000

    if insertion_rate > 0:
        assert len(aln) >= length
    else:
        assert len(aln) == length

    assert aln.num_seqs == len(tree.tips())


@pytest.mark.parametrize(
    "model",
    [
        *StandardDnaModel.iter_available_models(),
        *LieModel.iter_available_models(),
        *AaModel.iter_available_models(),
    ],
)
def test_rooted_tree(
    five_taxon_rooted_tree: PhyloNode,
    model: SubstitutionModel,
) -> None:
    check_simulate_alignment(five_taxon_rooted_tree, Model(model))


@pytest.mark.parametrize(
    "model",
    [
        *StandardDnaModel.iter_available_models(),
        *LieModel.iter_available_models(),
        *AaModel.iter_available_models(),
    ],
)
def test_unrooted_tree(
    four_taxon_unrooted_tree: PhyloNode,
    model: SubstitutionModel,
) -> None:
    check_simulate_alignment(four_taxon_unrooted_tree, Model(model))


@pytest.mark.parametrize("length", [None, 500, 1500])
def test_lengths(
    four_taxon_unrooted_tree: PhyloNode,
    length: int | None,
) -> None:
    check_simulate_alignment(
        four_taxon_unrooted_tree,
        Model(StandardDnaModel.GTR),
        length=length,
    )


@pytest.mark.parametrize("insertion_rate", [0.0, 0.03, 0.2])
@pytest.mark.parametrize("deletion_rate", [0.0, 0.04, 0.18])
def test_indel(
    four_taxon_unrooted_tree: PhyloNode,
    insertion_rate: float,
    deletion_rate: float,
) -> None:
    check_simulate_alignment(
        four_taxon_unrooted_tree,
        Model(StandardDnaModel.GTR),
        insertion_rate=insertion_rate,
        deletion_rate=deletion_rate,
    )
