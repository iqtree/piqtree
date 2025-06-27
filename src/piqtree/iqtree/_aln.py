"""Python wrapper for AliSim Simulation in the IQ-TREE library."""

import tempfile
from pathlib import Path

import cogent3
import cogent3.app.typing as c3_types
import yaml
from _piqtree import iq_simulate_alignment

from piqtree.iqtree._decorator import iqtree_func

iq_simulate_alignment = iqtree_func(iq_simulate_alignment, hide_files=True)


def simulate_alignment(
    trees: list[cogent3.PhyloNode],
    subst_model: str,
    seed: int,
    partition_info: list[str] | None = None,
    partition_type: str | None = None,
    seq_length: int | None = None,
    insertion_rate: float | None = None,
    deletion_rate: float | None = None,
    root_seq: str | None = None,
    num_threads: int | None = None,
    insertion_size_distribution: str | None = None,
    deletion_size_distribution: str | None = None,
) -> tuple[c3_types.AlignedSeqsType, str]:
    """Executes AliSim Simulation through IQ-TREE.

    Parameters
    ----------
    trees: list[cogent3.PhyloNode]
        A collection of trees.
    subst_model: str
        The substitution model name.
    seed: int
        The random seed.
    partition_info: list[str] | None, optional
        Partition information (by default None and will be set to []).
    partition_type: str | None, optional
        If provided, partition type must be 'equal', 'proportion', or 'unlinked' (by default None and will be set to "").
    seq_length: int | None, optional
        The length of sequences (by default None and will be set to 1000).
    insertion_rate: float | None, optional
        The insertion rate (by default None and will be set to 0.0).
    deletion_rate: float | None, optional
        The deletion rate (by default None and will be set to 0.0).
    root_seq: str | None, optional
        The root sequence (by default None and will be set to "").
    num_threads: int | None, optional
        The number of threads (by default None and will be set to 1).
    insertion_size_distribution: str | None, optional
        The insertion size distribution (by default None and will be set to "").
    deletion_size_distribution: str | None, optional
        The deletion size distribution (by default None and will be set to "").

    Returns
    -------
    tuple[c3_types.AlignedSeqsType, str]
        The simulated alignment and the content of the log file.

    """

    # Convert the trees to Newick strings
    newick_trees = [
        tree.get_newick(with_distances=True, semicolon=True) for tree in trees
    ]

    # Handle cases where some input variables are None
    if partition_info is None:
        partition_info = []
    if partition_type is None:
        partition_type = ""
    if seq_length is None:
        seq_length = 1000
    if insertion_rate is None:
        insertion_rate = 0.0
    if deletion_rate is None:
        deletion_rate = 0.0
    if root_seq is None:
        root_seq = ""
    if num_threads is None:
        num_threads = 1
    if insertion_size_distribution is None:
        insertion_size_distribution = ""
    if deletion_size_distribution is None:
        deletion_size_distribution = ""

    # Call the IQ-TREE function
    yaml_result = yaml.safe_load(
        iq_simulate_alignment(
            newick_trees,
            subst_model,
            seed,
            partition_info,
            partition_type,
            seq_length,
            insertion_rate,
            deletion_rate,
            root_seq,
            num_threads,
            insertion_size_distribution,
            deletion_size_distribution,
        ),
    )

    # Extract the simulated alignment and the content of the log file from the YAML result
    # cogent3.make_aligned_seqs dictionary
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
        tmp.write(yaml_result["alignment"])
        tmp.flush()
        aln = cogent3.load_aligned_seqs(
            tmp.name,
            format="phylip",
        )  # , new_type=True, moltype = 'text') #text for now
    Path(tmp.name).unlink()

    log_str = yaml_result["log"]

    return aln, log_str
