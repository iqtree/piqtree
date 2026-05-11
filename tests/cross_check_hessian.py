"""Cross-check piqtree fit_tree_hessian against iqtree3 --dating mcmctree.

For each (alignment, model) test case, both pipelines should produce the same
branch-length vector, gradient and Hessian matrix (modulo numerical noise)
because they go through the same underlying IQ-TREE computation.

Run with:
    python3 tests/cross_check_hessian.py
"""

from __future__ import annotations

import os
import pathlib
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent

# Make the freshly-built _piqtree extension importable without going through
# the full piqtree Python package (which needs Python 3.12+ and cogent3).
sys.path.insert(0, str(REPO / "src" / "piqtree" / "_libiqtree"))
import _piqtree  # noqa: E402

CLI = REPO / "iqtree3" / "build_bin" / "iqtree3"
RAND_SEED = 42


@dataclass(frozen=True)
class Case:
    name: str
    fasta: pathlib.Path
    taxa: list[str]
    tree_topology: str  # Newick — taxa names will be looked up in the alignment
    model: str


CASES = [
    Case(
        name="DNA / JC",
        fasta=REPO / "tests" / "data" / "example.fasta",
        taxa=["Human", "Chimpanzee", "Mouse", "Rhesus"],
        tree_topology="((Human,Chimpanzee),(Mouse,Rhesus));",
        model="JC",
    ),
    Case(
        name="DNA / HKY+G",
        fasta=REPO / "tests" / "data" / "example.fasta",
        taxa=["Human", "Chimpanzee", "Mouse", "Rhesus"],
        tree_topology="((Human,Chimpanzee),(Mouse,Rhesus));",
        model="HKY+G",
    ),
    Case(
        name="DNA / GTR+G",
        fasta=REPO / "tests" / "data" / "example.fasta",
        taxa=["Human", "Chimpanzee", "Mouse", "Rhesus"],
        tree_topology="((Human,Chimpanzee),(Mouse,Rhesus));",
        model="GTR+G",
    ),
    Case(
        name="AA / LG+G",
        fasta=REPO / "tests" / "data" / "protein.fasta",
        taxa=["Seq1", "Seq2", "Seq3", "Seq4"],
        tree_topology="((Seq1,Seq2),(Seq3,Seq4));",
        model="LG+G",
    ),
]


def read_fasta(path: pathlib.Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur: str | None = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].strip()
            seqs[cur] = ""
        elif cur is not None:
            seqs[cur] += line.strip()
    return seqs


def trim_gaps(rows: list[str], alphabet: str) -> list[str]:
    """Drop columns where any sequence has a gap or ambiguous residue.

    `alphabet` is "dna" or "aa"; ambiguity sets differ.
    """
    if alphabet == "dna":
        skip = set("-?NXnx")
    else:
        skip = set("-?XBZJUOxbzjuo")
    n = len(rows[0])
    keep = [i for i in range(n) if not any(r[i] in skip for r in rows)]
    return ["".join(r[i] for i in keep) for r in rows]


def parse_hessian_file(path: pathlib.Path) -> dict:
    """Parse an MCMCTree-style .mcmctree.hessian file into structured values.

    Layout (non-partitioned case):
        <nseq>
        <newick>
        <branch lengths>    (one row of floats)
        <gradient>          (one row of floats)
        Hessian
        <n rows of n floats>
    """
    text = path.read_text()
    blocks = [b.strip() for b in re.split(r"\n\s*\n", text) if b.strip()]

    nseq = int(blocks[0])
    tree = blocks[1]
    branch_lengths = np.fromstring(blocks[2], sep=" ")
    gradient = np.fromstring(blocks[3], sep=" ")

    hess_block = blocks[4]
    first_newline = hess_block.find("\n")
    label = hess_block[:first_newline].strip()
    assert label.lower().startswith("hessian"), label
    matrix_text = hess_block[first_newline + 1:]

    n = branch_lengths.size
    rows = np.fromstring(matrix_text, sep=" ")
    if rows.size != n * n:
        raise ValueError(
            f"hessian row count mismatch: got {rows.size} elements, "
            f"expected {n * n} ({n}x{n})",
        )
    hessian = rows.reshape(n, n)
    return {
        "nseq": nseq,
        "tree": tree,
        "branch_lengths": branch_lengths,
        "gradient": gradient,
        "hessian": hessian,
    }


def run_cli(workdir: pathlib.Path, names: list[str], seqs: list[str],
            tree_newick: str, model: str) -> dict:
    aln = workdir / "aln.fasta"
    with aln.open("w") as fh:
        for n, s in zip(names, seqs):
            fh.write(f">{n}\n{s}\n")
    tree = workdir / "start.nwk"
    tree.write_text(tree_newick + "\n")

    prefix = workdir / "cli"
    cmd = [
        str(CLI),
        "-s", str(aln),
        "-te", str(tree),
        "-m", model,
        "--dating", "mcmctree",
        "-nt", "1",
        "-seed", str(RAND_SEED),
        "--prefix", str(prefix),
        "-redo",
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return parse_hessian_file(workdir / "cli.mcmctree.hessian")


def run_piqtree(names: list[str], seqs: list[str], tree_newick: str,
                model: str, workdir: pathlib.Path) -> dict:
    # IQ-TREE writes some side files into cwd; isolate them so we don't pollute
    # the repo. Also silence its verbose stdout/stderr.
    old = pathlib.Path.cwd()
    os.chdir(workdir)
    try:
        devnull = os.open(os.devnull, os.O_WRONLY)
        saved_out = os.dup(1)
        saved_err = os.dup(2)
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        try:
            result = _piqtree.iq_fit_tree_hessian(
                names, seqs, model, tree_newick,
                False, RAND_SEED, 1, "",
            )
        finally:
            os.dup2(saved_out, 1)
            os.dup2(saved_err, 2)
            os.close(saved_out)
            os.close(saved_err)
            os.close(devnull)
    finally:
        os.chdir(old)
    return result


def compare(case: Case, cli: dict, piq: dict) -> bool:
    print(f"  branch_lengths CLI: {cli['branch_lengths']}")
    print(f"  branch_lengths piq: {piq['branch_lengths']}")
    print(f"  gradient       CLI: {cli['gradient']}")
    print(f"  gradient       piq: {piq['gradient']}")
    if cli["hessian"].size <= 49:  # only print small matrices in full
        with np.printoptions(precision=4, suppress=False, linewidth=140):
            print(f"  hessian        CLI:\n{cli['hessian']}")
            print(f"  hessian        piq:\n{piq['hessian']}")

    ok = True
    for key in ("branch_lengths", "gradient", "hessian"):
        a = np.asarray(cli[key], dtype=float)
        b = np.asarray(piq[key], dtype=float)
        if a.shape != b.shape:
            print(f"  FAIL [{key}]: shape mismatch CLI {a.shape} vs piqtree {b.shape}")
            ok = False
            continue
        diff = np.abs(a - b)
        max_abs = float(diff.max())
        denom = np.maximum(np.abs(a), 1e-12)
        max_rel = float((diff / denom).max())
        # CLI text output is ~6 sig figs, so use rtol=1e-3 and atol=1e-5.
        close = np.allclose(a, b, atol=1e-5, rtol=1e-3)
        print(f"  [{key}] max_abs={max_abs:.3e} max_rel={max_rel:.3e} -> "
              f"{'OK' if close else 'MISMATCH'}")
        if not close:
            ok = False
    return ok


def main() -> int:
    if not CLI.exists():
        print(f"ERROR: CLI binary not found at {CLI}")
        return 1

    overall_ok = True
    for case in CASES:
        print(f"\n=== {case.name} ===")
        full = read_fasta(case.fasta)
        try:
            selected = [full[t] for t in case.taxa]
        except KeyError as exc:
            print(f"  SKIP: taxon {exc!s} missing in {case.fasta.name}")
            overall_ok = False
            continue
        alphabet = "aa" if case.fasta.name == "protein.fasta" else "dna"
        cleaned = trim_gaps(selected, alphabet)
        if len({len(r) for r in cleaned}) != 1:
            print(f"  SKIP: ungapped alignment not square for {case.name}")
            overall_ok = False
            continue
        print(f"  alignment: {len(case.taxa)} taxa, {len(cleaned[0])} sites, "
              f"alphabet={alphabet}")

        slug = re.sub(r"[^A-Za-z0-9]+", "_", case.name).strip("_")
        with tempfile.TemporaryDirectory(prefix=f"hcmp_{slug}_") as td:
            workdir = pathlib.Path(td)
            (workdir / "cli").mkdir()
            (workdir / "piq").mkdir()
            cli = run_cli(workdir / "cli", case.taxa, cleaned, case.tree_topology, case.model)
            piq = run_piqtree(case.taxa, cleaned, case.tree_topology, case.model,
                              workdir / "piq")
        ok = compare(case, cli, piq)
        print(f"  -> {'PASS' if ok else 'FAIL'}")
        overall_ok &= ok

    print()
    print("=" * 60)
    print(f"OVERALL: {'PASS' if overall_ok else 'FAIL'}")
    return 0 if overall_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
