"""cogent3 apps for piqtree."""

from collections.abc import Iterable

import cogent3
import cogent3.app.typing as c3_types
from cogent3.app import composable
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
)
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
    ) -> None:
        self._model = model
        self._rand_seed = rand_seed
        self._bootstrap_reps = bootstrap_reps
        self._num_threads = num_threads

    def main(
        self,
        aln: c3_types.AlignedSeqsType,
    ) -> cogent3.PhyloNode | cogent3.app.typing.SerialisableType:
        tree = build_tree(
            aln,
            self._model,
            self._rand_seed,
            bootstrap_replicates=self._bootstrap_reps,
            num_threads=self._num_threads,
        )
        tree.source = getattr(aln, "source", None)
        return tree


@composable.define_app
class piq_fit_tree:
    @extend_docstring_from(fit_tree)
    def __init__(
        self,
        tree: cogent3.PhyloNode,
        model: Model | str,
        *,
        num_threads: int | None = None,
        bl_fixed: bool = False,
    ) -> None:
        self._tree = tree
        self._model = model
        self._num_threads = num_threads
        self._bl_fixed = bl_fixed

    def main(
        self,
        aln: c3_types.AlignedSeqsType,
    ) -> cogent3.PhyloNode | cogent3.app.typing.SerialisableType:
        tree = fit_tree(
            aln,
            self._tree,
            self._model,
            self._num_threads,
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
) -> cogent3.PhyloNode:
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
        aln: c3_types.AlignedSeqsType,
    ) -> c3_types.PairwiseDistanceType | cogent3.app.typing.SerialisableType:
        dists = jc_distances(
            aln,
            num_threads=self._num_threads,
        )
        dists.source = getattr(aln, "source", None)
        return dists


@composable.define_app
@extend_docstring_from(nj_tree)
def piq_nj_tree(
    dists: c3_types.PairwiseDistanceType,
    *,
    allow_negative: bool = False,
) -> cogent3.PhyloNode:
    tree = nj_tree(dists, allow_negative=allow_negative)
    tree.params |= {"provenance": "piqtree"}
    tree.source = getattr(dists, "source", None)
    return tree


@composable.define_app
@extend_docstring_from(model_finder)
def piq_model_finder(
    aln: c3_types.AlignedSeqsType,
) -> ModelFinderResult | c3_types.SerialisableType:
    return model_finder(aln)


@composable.define_app
@extend_docstring_from(consensus_tree)
def piq_consensus_tree(
    trees: Iterable[cogent3.PhyloNode],
    *,
    min_support: float = 0.5,
) -> cogent3.PhyloNode:
    return consensus_tree(trees, min_support=min_support)


_ALL_APP_NAMES = [
    "piq_build_tree",
    "piq_fit_tree",
    "piq_random_tree",
    "piq_jc_distances",
    "piq_nj_tree",
    "piq_model_finder",
    "piq_consensus_tree",
]
