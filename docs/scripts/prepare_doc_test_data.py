from pathlib import Path

from cogent3 import load_aligned_seqs

PROJECT_ROOT = Path(__file__)
while PROJECT_ROOT.name != "docs":
    PROJECT_ROOT = PROJECT_ROOT.parent
PROJECT_ROOT = PROJECT_ROOT.parent


def prepare_my_alignment() -> None:
    aln = load_aligned_seqs(PROJECT_ROOT / "tests/data/example.fasta")
    aln = aln.take_seqs(["Human", "Chimpanzee", "Rhesus", "Mouse"])
    aln.write("my_alignment.fasta")


def prepare_protein_alignment() -> None:
    aln = load_aligned_seqs(PROJECT_ROOT / "tests/data/protein.fasta")
    aln = aln.take_seqs(sorted(aln.names)[:4])
    aln.write("my_protein.fasta")


if __name__ == "__main__":
    prepare_my_alignment()
    prepare_protein_alignment()
