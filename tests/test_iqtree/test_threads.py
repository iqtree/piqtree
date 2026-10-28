import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest
from cogent3 import make_aligned_seqs
from cogent3.core.alignment import Alignment

from piqtree import build_tree, jc_distances, nj_tree


def _related_alignment(num_seqs: int, length: int) -> Alignment:
    rng = np.random.default_rng(0)
    bases = np.frombuffer(b"ACGT", dtype=np.uint8)
    root = rng.integers(4, size=length)
    seqs = {}
    for i in range(num_seqs):
        seq = root.copy()
        mutate = rng.random(length) < 0.05
        seq[mutate] = rng.integers(4, size=int(mutate.sum()))
        seqs[f"s{i}"] = bases[seq].tobytes().decode()
    return make_aligned_seqs(seqs, moltype="dna")


def test_other_threads_run_during_iqtree_call() -> None:
    # The GIL is released while IQ-TREE computes, so this thread keeps running.
    # If it were held, this thread would stall for the whole IQ-TREE call.
    aln = _related_alignment(num_seqs=20, length=500_000)
    started = threading.Event()
    done = threading.Event()

    def compute() -> None:
        started.set()
        jc_distances(aln, num_threads=1)
        done.set()

    worker = threading.Thread(target=compute)
    worker.start()
    started.wait()

    start = last = time.perf_counter()
    max_gap = 0.0
    while not done.is_set():
        now = time.perf_counter()
        max_gap = max(max_gap, now - last)
        last = now
    duration = last - start
    worker.join()

    assert max_gap < duration / 3


def test_concurrent_calls(
    five_otu: Alignment,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Calls from several threads give the same results as sequential calls, and
    # IQ-TREE's output files never end up in the current directory.
    monkeypatch.chdir(tmp_path)

    def run(i: int) -> str:
        if i % 2:
            tree = build_tree(five_otu, "JC", rand_seed=1)
        else:
            tree = nj_tree(jc_distances(five_otu))
        return tree.get_newick(with_distances=True)

    expected = [run(i) for i in range(2)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        actual = list(pool.map(run, range(16)))

    assert actual == expected * 8
    assert list(tmp_path.iterdir()) == []
