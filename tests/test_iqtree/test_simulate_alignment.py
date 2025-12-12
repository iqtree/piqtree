import re
from typing import cast

import pytest
from cogent3 import make_tree
from cogent3.core.alignment import Alignment
from cogent3.core.tree import PhyloNode

from piqtree import Model, simulate_alignment
from piqtree.distribution import (
    IndelDistribution,
    IndelGeometric,
    IndelLavalette,
    IndelNegativeBinomial,
    IndelZipfian,
)
from piqtree.iqtree._alignment import UNSUPPORTED_MODELS
from piqtree.model import (
    AaModel,
    LieModel,
    StandardDnaModel,
    SubstitutionModel,
)


def check_simulate_alignment(
    tree: PhyloNode,
    model: Model | str,
    length: int | None = None,
    insertion_rate: float = 0.0,
    deletion_rate: float = 0.0,
    insertion_size_distribution: IndelDistribution | str = "POW{1.7/100}",
    deletion_size_distribution: IndelDistribution | str = "POW{1.7/100}",
    seed: int | None = None,
) -> Alignment:
    if length is None:
        aln = simulate_alignment(
            tree,
            model,
            rand_seed=seed,
            insertion_rate=insertion_rate,
            deletion_rate=deletion_rate,
            insertion_size_distribution=insertion_size_distribution,
            deletion_size_distribution=deletion_size_distribution,
        )
    else:
        aln = simulate_alignment(
            tree,
            model,
            length,
            rand_seed=seed,
            insertion_rate=insertion_rate,
            deletion_rate=deletion_rate,
            insertion_size_distribution=insertion_size_distribution,
            deletion_size_distribution=deletion_size_distribution,
        )
    if length is None:
        length = 1000

    if insertion_rate > 0:
        assert len(aln) >= length
    else:
        assert len(aln) == length

    assert aln.num_seqs == len(tree.tips())
    assert sorted(aln.names) == sorted(tree.get_tip_names())

    return aln


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
    if model.base_model in UNSUPPORTED_MODELS:
        with pytest.raises(
            ValueError,
            match=re.escape(
                f"Lie Model {cast('LieModel', model.base_model).value} is unsupported.",
            ),
        ):
            check_simulate_alignment(five_taxon_rooted_tree, Model(model))
    else:
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
    if model.base_model in UNSUPPORTED_MODELS:
        with pytest.raises(
            ValueError,
            match=re.escape(
                f"Lie Model {cast('LieModel', model.base_model).value} is unsupported.",
            ),
        ):
            check_simulate_alignment(four_taxon_unrooted_tree, Model(model))
    else:
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
def test_indel_rates(
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


@pytest.mark.parametrize(
    "indel_size_distributions",
    [
        (IndelGeometric(5), IndelGeometric(4)),
        (IndelNegativeBinomial(5, 20), IndelZipfian(1.5, 10)),
        (IndelLavalette(1.5, 10), IndelLavalette(1.5, 10)),
    ],
)
def test_indel_sizes(
    four_taxon_unrooted_tree: PhyloNode,
    indel_size_distributions: tuple[IndelDistribution | str, IndelDistribution | str],
) -> None:
    insertion_size_distribution = indel_size_distributions[0]
    deletion_size_distribution = indel_size_distributions[1]
    aln = check_simulate_alignment(
        four_taxon_unrooted_tree,
        "GTR{1.0,2.0,1.5,3.7,2.8}+F{0.1,0.2,0.3,0.4}",
        insertion_rate=0.1,
        deletion_rate=0.05,
        insertion_size_distribution=insertion_size_distribution,
        deletion_size_distribution=deletion_size_distribution,
        seed=1,
    )

    aln_indel_str = check_simulate_alignment(
        four_taxon_unrooted_tree,
        "GTR{1.0,2.0,1.5,3.7,2.8}+F{0.1,0.2,0.3,0.4}",
        insertion_rate=0.1,
        deletion_rate=0.05,
        insertion_size_distribution=str(insertion_size_distribution),
        deletion_size_distribution=str(deletion_size_distribution),
        seed=1,
    )

    assert aln.names == aln_indel_str.names

    for name in aln.names:
        assert aln.get_seq(name) == aln_indel_str.get_seq(name)


def test_root_seq() -> None:
    tree = make_tree("((a:0.8,(b:1.2,c:0)):0.0,d:0.5,e:0.0)")
    root_seq = "GGGGCCCCAAAATTTT" * 10

    aln = simulate_alignment(tree, "JC", len(root_seq), root_seq=root_seq)

    assert str(aln.get_seq("c")) == root_seq
    assert str(aln.get_seq("e")) == root_seq
