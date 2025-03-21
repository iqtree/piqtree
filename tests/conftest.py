import pathlib

import pytest
from cogent3 import Alignment, load_aligned_seqs


@pytest.fixture(scope="session")
def DATA_DIR() -> pathlib.Path:
    return pathlib.Path(__file__).parent / "data"


@pytest.fixture
def three_otu(DATA_DIR: pathlib.Path) -> Alignment:
    aln = load_aligned_seqs(DATA_DIR / "example.fasta", moltype="dna", new_type=True)
    aln = aln.take_seqs(["Human", "Rhesus", "Mouse"])
    return aln.omit_gap_pos(allowed_gap_frac=0)


@pytest.fixture
def four_otu(DATA_DIR: pathlib.Path) -> Alignment:
    aln = load_aligned_seqs(DATA_DIR / "example.fasta", moltype="dna", new_type=True)
    aln = aln.take_seqs(["Human", "Chimpanzee", "Rhesus", "Mouse"])
    return aln.omit_gap_pos(allowed_gap_frac=0)


@pytest.fixture
def five_otu(DATA_DIR: pathlib.Path) -> Alignment:
    aln = load_aligned_seqs(DATA_DIR / "example.fasta", moltype="dna", new_type=True)
    aln = aln.take_seqs(["Human", "Chimpanzee", "Rhesus", "Manatee", "Dugong"])
    return aln.omit_gap_pos(allowed_gap_frac=0)


@pytest.fixture
def all_otu(DATA_DIR: pathlib.Path) -> Alignment:
    aln = load_aligned_seqs(DATA_DIR / "example.fasta", moltype="dna", new_type=True)
    return aln.omit_gap_pos(allowed_gap_frac=0)
