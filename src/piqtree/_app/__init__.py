"""cogent3 apps for piqtree."""

from collections.abc import Iterable

from cogent3.app import composable
from cogent3.core.alignment import Alignment
from cogent3.core.tree import PhyloNode
from cogent3.evolve.fast_distance import DistanceMatrix
from cogent3.util.misc import extend_docstring_from

from piqtree import (
    TreeGenMode,
    build_tree,
    consensus_tree,
    fit_tree,
    jc_distances,
    model_finder,
    nj_tree,
    random_tree,
    simulate_alignment,
)
from piqtree.distribution import IndelDistribution
from piqtree.iqtree import ModelFinderResult
from piqtree.model import Model


@composable.define_app
class piq_build_tree:
    @extend_docstring_from(build_tree)
    def __init__(
        self,
        model: Model | str,
        *,
        rand_seed: int | None = None,
        bootstrap_reps: int | None = None,
        num_threads: int | None = None,
        other_options: str = "",
    ) -> None:
        self._model = model
        self._rand_seed = rand_seed
        self._bootstrap_reps = bootstrap_reps
        self._num_threads = num_threads
        self._other_options = other_options

    def main(
        self,
        aln: Alignment,
    ) -> PhyloNode:
        tree = build_tree(
            aln,
            self._model,
            self._rand_seed,
            bootstrap_replicates=self._bootstrap_reps,
            num_threads=self._num_threads,
            other_options=self._other_options,
        )
        tree.source = getattr(aln, "source", None)
        return tree


@composable.define_app
class piq_fit_tree:
    @extend_docstring_from(fit_tree)
    def __init__(
        self,
        tree: PhyloNode,
        model: Model | str,
        *,
        num_threads: int | None = None,
        other_options: str = "",
        bl_fixed: bool = False,
    ) -> None:
        self._tree = tree
        self._model = model
        self._num_threads = num_threads
        self._bl_fixed = bl_fixed
        self._other_options = other_options

    def main(
        self,
        aln: Alignment,
    ) -> PhyloNode:
        tree = fit_tree(
            aln,
            self._tree,
            self._model,
            self._num_threads,
            other_options=self._other_options,
            bl_fixed=self._bl_fixed,
        )
        tree.source = getattr(aln, "source", None)
        return tree


@composable.define_app
@extend_docstring_from(random_tree)
def piq_random_tree(
    num_taxa: int,
    tree_mode: TreeGenMode,
    rand_seed: int | None = None,
) -> PhyloNode:
    return random_tree(num_taxa, tree_mode, rand_seed)


@composable.define_app
class piq_jc_distances:
    @extend_docstring_from(jc_distances)
    def __init__(
        self,
        num_threads: int | None = None,
    ) -> None:
        self._num_threads = num_threads

    def main(
        self,
        aln: Alignment,
    ) -> DistanceMatrix:
        dists = jc_distances(
            aln,
            num_threads=self._num_threads,
        )
        dists.source = getattr(aln, "source", None)
        return dists


@composable.define_app
@extend_docstring_from(nj_tree)
def piq_nj_tree(
    dists: DistanceMatrix,
    *,
    allow_negative: bool = False,
) -> PhyloNode:
    tree = nj_tree(dists, allow_negative=allow_negative)
    tree.params |= {"provenance": "piqtree"}
    tree.source = getattr(dists, "source", None)
    return tree


@composable.define_app
@extend_docstring_from(model_finder)
def piq_model_finder(
    aln: Alignment,
) -> ModelFinderResult:
    return model_finder(aln)


@composable.define_app
@extend_docstring_from(consensus_tree)
def piq_consensus_tree(
    trees: Iterable[PhyloNode],
    *,
    min_support: float = 0.5,
) -> PhyloNode:
    return consensus_tree(trees, min_support=min_support)


@composable.define_app
@extend_docstring_from(simulate_alignment)
def piq_simulate_alignment(
    tree: PhyloNode,
    model: Model | str,
    length: int = 1000,
    rand_seed: int | None = None,
    insertion_rate: float = 0.0,
    deletion_rate: float = 0.0,
    insertion_size_distribution: IndelDistribution | str = "POW{1.7/100}",
    deletion_size_distribution: IndelDistribution | str = "POW{1.7/100}",
    root_seq: str | None = None,
    num_threads: int | None = None,
) -> Alignment:
    return simulate_alignment(
        tree,
        model,
        length=length,
        rand_seed=rand_seed,
        insertion_rate=insertion_rate,
        deletion_rate=deletion_rate,
        insertion_size_distribution=insertion_size_distribution,
        deletion_size_distribution=deletion_size_distribution,
        root_seq=root_seq,
        num_threads=num_threads,
    )


_ALL_APP_NAMES = [
    "piq_build_tree",
    "piq_fit_tree",
    "piq_random_tree",
    "piq_jc_distances",
    "piq_nj_tree",
    "piq_model_finder",
    "piq_consensus_tree",
    "piq_simulate_alignment",
]
