import os
import tempfile
from pathlib import Path

import nox

_py_versions = range(13, 16)
_python_sessions = [f"3.{v}" for v in _py_versions]

nox.options.default_venv_backend = "uv"
# Always use uv-managed Python builds, rather than whatever interpreter happens
# to be on the machine, so sessions run on the same builds locally and in CI.
nox.options.download_python = "always"


def _uv_sync(session: nox.Session, *args: str) -> None:
    """Install piqtree and locked dependencies from uv.lock into the session venv."""
    session.run_install(
        "uv",
        "sync",
        "--locked",
        "--no-default-groups",
        f"--python={session.virtualenv.location}",
        *args,
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )


@nox.session(python=_python_sessions)
def test(session: nox.Session) -> None:
    posargs = list(session.posargs)
    env = os.environ.copy()

    _uv_sync(session, "--group", "test")
    session.run("pytest", *posargs, env=env)


@nox.session(python=_python_sessions)
def type_check(session: nox.Session) -> None:
    posargs = list(session.posargs)
    env = os.environ.copy()

    _uv_sync(session, "--no-editable", "--group", "typing")
    _check_stub(session)
    session.run("mypy", "src", "tests", *posargs, env=env)


_STUB = Path("src/piqtree/_piqtree.pyi")
_STUBGEN = (
    "python",
    "-m",
    "nanobind.stubgen",
    "-q",
    "--exclude-values",
    "-m",
    "piqtree._piqtree",
)


def _check_stub(session: nox.Session) -> None:
    """Check the committed stub for the C++ extension is up to date."""
    with tempfile.TemporaryDirectory() as tmpdir:
        generated = Path(tmpdir) / _STUB.name
        session.run(*_STUBGEN, "-o", str(generated))
        if generated.read_text() != _STUB.read_text():
            session.error(
                f"{_STUB} is out of date. Regenerate it by running:\n"
                f"    uv run {' '.join(_STUBGEN)} -o {_STUB}",
            )


@nox.session(python=_python_sessions)
def ruff(session: nox.Session) -> None:
    posargs = list(session.posargs)
    env = os.environ.copy()

    _uv_sync(session, "--no-editable", "--group", "lint")
    session.run("ruff", "check", *posargs, env=env)


@nox.session(python=_python_sessions)
def test_docs(session: nox.Session) -> None:
    doc_md_files = [path.resolve() for path in Path("docs").rglob("*.md")]
    doc_py_files = [path.resolve() for path in Path("docs").rglob("*.py")]

    doctest_setup_path = Path("docs/scripts/prepare_doc_test_data.py").resolve()

    posargs: list[str | Path] = list(session.posargs)
    posargs.extend(["--markdown-docs", *doc_md_files])
    env = os.environ.copy()

    _uv_sync(session, "--no-editable", "--group", "test")

    with tempfile.TemporaryDirectory() as tmpdir:
        session.chdir(tmpdir)

        session.run("python", doctest_setup_path, env=env)
        session.run("pytest", *posargs, env=env)

        session.install("cogent3[extra]", "diverse-seq")
        for py_file in doc_py_files:
            session.run("python", py_file, env=env, silent=True)
