import pathlib

import pytest
from cogent3 import PhyloNode, load_aligned_seqs, make_tree
from cogent3.core.alignment import Alignment


@pytest.fixture(scope="session")
def DATA_DIR() -> pathlib.Path:
    return pathlib.Path(__file__).parent / "data"


@pytest.fixture
def three_otu(DATA_DIR: pathlib.Path) -> Alignment:
    aln = load_aligned_seqs(DATA_DIR / "example.fasta", moltype="dna")
    aln = aln.take_seqs(["Human", "Rhesus", "Mouse"])
    return aln.omit_gap_pos(allowed_gap_frac=0)


@pytest.fixture
def four_otu(DATA_DIR: pathlib.Path) -> Alignment:
    aln = load_aligned_seqs(
        DATA_DIR / "example.fasta",
        moltype="dna",
    )
    aln = aln.take_seqs(["Human", "Chimpanzee", "HumpbackW", "SpermWhale"])
    aln = aln.omit_gap_pos(allowed_gap_frac=0)
    return aln[::15]


@pytest.fixture
def five_otu(DATA_DIR: pathlib.Path) -> Alignment:
    aln = load_aligned_seqs(
        DATA_DIR / "example.fasta",
        moltype="dna",
    )
    aln = aln.take_seqs(["Human", "Chimpanzee", "Rhesus", "Manatee", "Dugong"])
    return aln.omit_gap_pos(allowed_gap_frac=0)


@pytest.fixture
def all_otu(DATA_DIR: pathlib.Path) -> Alignment:
    aln = load_aligned_seqs(
        DATA_DIR / "example.fasta",
        moltype="dna",
    )
    return aln.omit_gap_pos(allowed_gap_frac=0)


@pytest.fixture
def protein_four_otu(DATA_DIR: pathlib.Path) -> Alignment:
    aln = load_aligned_seqs(DATA_DIR / "protein.fasta", moltype="protein")
    aln = aln.take_seqs(sorted(aln.names)[:4])
    aln = aln.omit_gap_pos(allowed_gap_frac=0)
    return aln[::15]


@pytest.fixture
def five_trees() -> list[PhyloNode]:
    tree1 = make_tree("(a,(b,(c,(d,(e,f)))))")
    tree2 = make_tree("(a,(b,(c,(d,(e,f)))))")
    tree3 = make_tree("((a,b),(c,(d,(e,f))))")
    tree4 = make_tree("(((a,b),c),(d,(e,f)))")
    tree5 = make_tree("((((a,b),c),d),(e,f))")

    return [tree1, tree2, tree3, tree4, tree5]


@pytest.fixture
def four_taxon_unrooted_tree() -> PhyloNode:  # Unrooted
    return make_tree("((a:0.1,b:0.2):0.05,c:0.3,d:0.1);")


@pytest.fixture
def five_taxon_rooted_tree() -> PhyloNode:  # Rooted
    return make_tree("(((a:0.1,b:0.2):0.05,(c:0.3,d:0.1):0.2):0.05,e:0.4);")
