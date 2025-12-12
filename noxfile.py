import os
import tempfile
from pathlib import Path

import nox

_py_versions = range(12, 15)
_python_sessions = [f"3.{v}" for v in _py_versions]


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
    doc_md_files = [path.resolve() for path in Path("docs").rglob("*.md")]
    doc_py_files = [path.resolve() for path in Path("docs").rglob("*.py")]

    doctest_setup_path = Path("docs/scripts/prepare_doc_test_data.py").resolve()

    posargs: list[str | Path] = list(session.posargs)
    posargs.extend(["--markdown-docs", *doc_md_files])
    env = os.environ.copy()

    install_spec = ".[test]"
    session.install(install_spec)

    with tempfile.TemporaryDirectory() as tmpdir:
        session.chdir(tmpdir)

        session.run("python", doctest_setup_path, env=env)
        session.run("pytest", *posargs, env=env)

        session.install("cogent3[extra]", "diverse-seq")
        for py_file in doc_py_files:
            session.run("python", py_file, env=env, silent=True)
