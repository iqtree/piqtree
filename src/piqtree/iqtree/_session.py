"""Persistent IQ-TREE session for repeated likelihood / Hessian queries."""

from __future__ import annotations

import os
import pathlib
import tempfile
from dataclasses import dataclass
from typing import Any, cast

import numpy as np
from _piqtree import (
    iq_hessian_session_create,
    iq_hessian_session_eval,
    iq_hessian_session_score,
)
from cogent3 import make_tree
from cogent3.core.alignment import Alignment
from cogent3.core.tree import PhyloNode

from piqtree.iqtree._tree import _rename_iq_tree
from piqtree.model import Model, make_model
from piqtree.util import get_newick, process_rand_seed_nonzero


@dataclass(frozen=True)
class HessianAt:
    log_likelihood: float
    gradient: np.ndarray
    hessian: np.ndarray


def _tree_has_all_branch_lengths(tree: PhyloNode) -> bool:
    return all(node.length is not None for node in tree.preorder(include_self=False))


class IqtreeSession:
    """Persistent IQ-TREE session for repeated Hessian / likelihood queries."""

    def __init__(
        self,
        aln: Alignment,
        tree: PhyloNode,
        model: Model | str = "GTR+G",
        rand_seed: int | None = None,
        num_threads: int | None = None,
        *,
        bl_fixed: bool = False,
    ) -> None:
        self._aln = aln
        self._names = aln.names
        self._seqs = [str(s) for s in aln.iter_seqs(self._names)]
        self._tree_obj = tree
        self._model_obj: Model = model if isinstance(model, Model) else make_model(model)
        self._bl_fixed = bl_fixed
        self._rand_seed = process_rand_seed_nonzero(rand_seed)
        self._num_threads = 1 if num_threads is None else num_threads

        self._workdir = tempfile.mkdtemp(prefix="piqtree_session_")
        self._capsule: Any = None
        self._branch_lengths: np.ndarray = np.zeros(0)
        self._open()

    def _open(self) -> None:
        newick = get_newick(self._tree_obj)
        cwd = pathlib.Path.cwd()
        os.chdir(self._workdir)
        try:
            self._capsule = iq_hessian_session_create(
                self._names, self._seqs, str(self._model_obj), newick,
                self._bl_fixed, self._rand_seed, self._num_threads, "",
            )
            r = iq_hessian_session_eval(self._capsule, np.array([], dtype=np.float64))
            self._branch_lengths = cast("np.ndarray", r["branch_lengths"])
            self._fitted_newick = cast("str", r["tree"])
        finally:
            os.chdir(cwd)

    def _reopen(self) -> None:
        self._capsule = None
        self._open()

    def close(self) -> None:
        self._capsule = None
        wd = getattr(self, "_workdir", None)
        if wd:
            import shutil
            shutil.rmtree(wd, ignore_errors=True)
            self._workdir = ""

    def __enter__(self) -> "IqtreeSession":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:  # noqa: BLE001
            pass

    @property
    def branch_lengths(self) -> np.ndarray:
        return self._branch_lengths

    @property
    def tree(self) -> PhyloNode:
        t = make_tree(self._fitted_newick)
        _rename_iq_tree(t, self._names)
        t.name_unnamed_nodes()
        t.params["model"] = str(self._model_obj)
        return t

    @property
    def model(self) -> str:
        return str(self._model_obj)

    def update_branchlens(self, branch_lengths: np.ndarray) -> None:
        bls = np.ascontiguousarray(branch_lengths, dtype=np.float64)
        if bls.size != self._branch_lengths.size:
            msg = (
                f"branch_lengths size mismatch: got {bls.size}, "
                f"expected {self._branch_lengths.size}"
            )
            raise ValueError(msg)
        cwd = pathlib.Path.cwd()
        os.chdir(self._workdir)
        try:
            iq_hessian_session_score(self._capsule, bls)
        finally:
            os.chdir(cwd)
        self._branch_lengths = bls

    def update_tree(self, tree: PhyloNode) -> None:
        self._tree_obj = tree
        self._reopen()
        if _tree_has_all_branch_lengths(tree):
            try:
                bls = np.array(
                    [n.length for n in tree.preorder(include_self=False)],
                    dtype=np.float64,
                )
                if bls.size == self._branch_lengths.size:
                    self.update_branchlens(bls)
            except (TypeError, ValueError):
                pass

    def update_model(self, model: Model | str) -> None:
        self._model_obj = model if isinstance(model, Model) else make_model(model)
        self._reopen()

    def score_at(
        self,
        *,
        branch_lengths: np.ndarray | None = None,
        tree: PhyloNode | None = None,
        model: Model | str | None = None,
    ) -> float:
        self._apply_kwargs(branch_lengths=branch_lengths, tree=tree, model=model)
        cwd = pathlib.Path.cwd()
        os.chdir(self._workdir)
        try:
            return float(iq_hessian_session_score(
                self._capsule, np.array([], dtype=np.float64),
            ))
        finally:
            os.chdir(cwd)

    def hessian_at(
        self,
        *,
        branch_lengths: np.ndarray | None = None,
        tree: PhyloNode | None = None,
        model: Model | str | None = None,
    ) -> HessianAt:
        self._apply_kwargs(branch_lengths=branch_lengths, tree=tree, model=model)
        cwd = pathlib.Path.cwd()
        os.chdir(self._workdir)
        try:
            r = iq_hessian_session_eval(self._capsule, np.array([], dtype=np.float64))
        finally:
            os.chdir(cwd)
        self._branch_lengths = cast("np.ndarray", r["branch_lengths"])
        self._fitted_newick = cast("str", r["tree"])
        return HessianAt(
            log_likelihood=float(cast("float", r["log_likelihood"])),
            gradient=cast("np.ndarray", r["gradient"]),
            hessian=cast("np.ndarray", r["hessian"]),
        )

    def _apply_kwargs(
        self,
        *,
        branch_lengths: np.ndarray | None,
        tree: PhyloNode | None,
        model: Model | str | None,
    ) -> None:
        if model is not None:
            self.update_model(model)
        if tree is not None:
            self.update_tree(tree)
        if branch_lengths is not None:
            self.update_branchlens(branch_lengths)


def iqtree_session_create(
    aln: Alignment,
    tree: PhyloNode,
    model: Model | str = "GTR+G",
    rand_seed: int | None = None,
    num_threads: int | None = None,
    *,
    bl_fixed: bool = False,
) -> IqtreeSession:
    return IqtreeSession(
        aln, tree, model, rand_seed=rand_seed, num_threads=num_threads,
        bl_fixed=bl_fixed,
    )
