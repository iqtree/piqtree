from cogent3 import load_aligned_seqs


def prepare_my_alignment() -> None:
    aln = load_aligned_seqs("tests/data/example.fasta")
    aln = aln.take_seqs(["Human", "Chimpanzee", "Rhesus", "Mouse"])
    aln.write("my_alignment.fasta")


def prepare_protein_alignment() -> None:
    aln = load_aligned_seqs("tests/data/protein.fasta")
    aln = aln.take_seqs(sorted(aln.names)[:4])
    aln.write("my_protein.fasta")


if __name__ == "__main__":
    prepare_my_alignment()
    prepare_protein_alignment()
