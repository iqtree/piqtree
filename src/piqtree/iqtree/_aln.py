"""Python wrapper for AliSim Simulation in the IQ-TREE library."""

import tempfile
from pathlib import Path
from typing import Literal

import cogent3
import cogent3.app.typing as c3_types
import yaml
from _piqtree import iq_simulate_alignment

from piqtree.distribution import IndelDistribution
from piqtree.iqtree._decorator import iqtree_func
from piqtree.model import Model

iq_simulate_alignment = iqtree_func(iq_simulate_alignment, hide_files=True)


def simulate_alignment(
    trees: list[cogent3.PhyloNode],
    model: Model | str,
    length: int = 1000,
    rand_seed: int | None = None,
    insertion_rate: float = 0.0,
    deletion_rate: float = 0.0,
    insertion_size_distribution: IndelDistribution | str = "POW{1.7/100}",
    deletion_size_distribution: IndelDistribution | str = "POW{1.7/100}",
    root_seq: str | None = None,
    partition_info: list[str] | None = None,
    partition_type: Literal["equal", "proportion", "unlinked"] | None = None,
    num_threads: int | None = None,
) -> tuple[c3_types.AlignedSeqsType, str]:
    """Executes AliSim Simulation through IQ-TREE.

    Parameters
    ----------
    trees: list[cogent3.PhyloNode]
        A collection of trees.
    subst_model: str
        The substitution model name.
    length: int
        The length of sequences (by default 1000).
    rand_seed : int | None, optional
        The random seed - 0 or None means no seed, by default None.
    insertion_rate: float | None, optional
        The insertion rate (by default 0.0).
    deletion_rate: float | None, optional
        The deletion rate (by default 0.0).
    insertion_size_distribution: str | None, optional
        The insertion size distribution (by default the Zipfian
        distribution with a=1.7 and maximum size 100).
    deletion_size_distribution: str | None, optional
        The deletion size distribution (by default the Zipfian
        distribution with a=1.7 and maximum size 100).
    root_seq: str | None, optional
        The root sequence (by default None).
    partition_info: list[str] | None, optional
        Partition information (by default None).
    partition_type: Literal["equal", "proportion", "unlinked"] | None, optional
        If provided, partition type must be "equal", "proportion", or
        "unlinked" (by default None).
    num_threads: int | None, optional
        Number of threads for IQ-TREE to use, by default None (single-threaded).

    Returns
    -------
    c3_types.AlignedSeqsType
        The simulated alignment.

    """
    if rand_seed is None:
        rand_seed = 0

    if root_seq is None:
        root_seq = ""

    if partition_info is None:
        partition_info = []
    if partition_type is None:
        partition_type = ""

    if num_threads is None:
        num_threads = 1

    model_str = str(model) if isinstance(model, Model) else model

    # Convert the trees to Newick strings
    newick_trees = [str(tree) for tree in trees]

    # Call the IQ-TREE function
    yaml_result = yaml.safe_load(
        iq_simulate_alignment(
            newick_trees,
            model_str,
            rand_seed,
            partition_info,
            partition_type,
            length,
            insertion_rate,
            deletion_rate,
            root_seq,
            num_threads,
            str(insertion_size_distribution),
            str(deletion_size_distribution),
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
            new_type=True,
        )  # moltype = 'dna' or 'protein')
    Path(tmp.name).unlink()

    return aln
