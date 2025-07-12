import pytest
import cogent3
from cogent3.core.alignment import Alignment

import piqtree

@pytest.fixture
def jc_model() -> str:
    return "JC"

@pytest.mark.parametrize(
    "tree",
    [
        cogent3.make_tree("((A:0.1,B:0.2):0.1,(C:0.3,D:0.4):0.2,E:0.5);"),
        cogent3.make_tree("(((F:0.2,E:0.4):0.3,(C:0.1,B:0.5):0.2):0.4,(A:0.3,D:0.2):0.5);"),
    ],
)
@pytest.mark.parametrize("seed", [1, 3])
@pytest.mark.parametrize("seq_length", [None, 2000])
@pytest.mark.parametrize("insertion_rate", [None, 0.5])
@pytest.mark.parametrize("deletion_rate", [None, 0.5])
@pytest.mark.parametrize("num_threads", [None, 5])
def test_simulate_alignment_single_tree_jc(
    tree: cogent3.PhyloNode,
    jc_model: str,
    seed: int,
    seq_length: int | None,
    insertion_rate: float | None,
    deletion_rate: float | None,
    num_threads: int | None,
) -> None:
    res = piqtree.simulate_alignment(
        trees = [tree],
        subst_model = jc_model,
        seed = seed,
        seq_length = seq_length,
        insertion_rate = insertion_rate,
        deletion_rate = deletion_rate,
        num_threads = num_threads,
    )
        
    # Checks if simulate_alignment is returning a tuple with two values
    assert isinstance(res, tuple) and len(res) == 2
    
    aln = res[0]
    
    # Checks if the number of sequences in the output alignment is equal to the number of taxa in the input tree
    assert aln.num_seqs == len(tree.tips())
    
    # If no insertion or deletion rates are given, checks if every sequence length in the output alignment is equal to the input sequence length
    if insertion_rate is None and deletion_rate is None:
        unique_seq_lengths = list(set(aln.get_lengths().to_dict().values()))
        actual_seq_length = seq_length if seq_length is not None else 1000
        assert len(unique_seq_lengths) == 1 and unique_seq_lengths[0] == actual_seq_length
