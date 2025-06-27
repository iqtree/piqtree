import os
import shutil
from pathlib import Path

import nox

_py_versions = range(11, 14)
_python_sessions = [f"3.{v}" + (".5" if v == 13 else "") for v in _py_versions]


@nox.session(python=_python_sessions)
def test(session: nox.Session) -> None:
    posargs = list(session.posargs)
    env = os.environ.copy()

    install_spec = "-e.[test]"
    session.install(install_spec)
    session.run("pytest", *posargs, env=env)


@nox.session(python=_python_sessions)
def type_check(session: nox.Session) -> None:
    posargs = list(session.posargs)
    env = os.environ.copy()

    install_spec = ".[typing]"
    session.install(install_spec)
    session.run("mypy", "src", "tests", *posargs, env=env)


@nox.session(python=_python_sessions)
def ruff(session: nox.Session) -> None:
    posargs = list(session.posargs)
    env = os.environ.copy()

    install_spec = ".[lint]"
    session.install(install_spec)
    session.run("ruff", "check", *posargs, env=env)


@nox.session(python=_python_sessions)
def test_docs(session: nox.Session) -> None:
    doc_md_files = Path("docs").rglob("*.md")
    doc_py_files = Path("docs").rglob("*.py")

    posargs = list(session.posargs)
    posargs.extend(["--markdown-docs", *doc_md_files])
    env = os.environ.copy()

    install_spec = ".[test]"
    session.install(install_spec)

    try:
        session.run("python", "docs/scripts/prepare_doc_test_data.py", env=env)
        session.run("pytest", *posargs, env=env)

        session.install("cogent3[extra]", "diverse-seq")
        for py_file in doc_py_files:
            session.run("python", py_file, env=env, silent=True)
    finally:
        Path("my_alignment.fasta").unlink(missing_ok=True)
        Path("my_protein.fasta").unlink(missing_ok=True)

        shutil.rmtree("data")
