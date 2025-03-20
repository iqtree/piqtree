import os

import nox

_py_versions = range(11, 14)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    posargs = list(session.posargs)
    env = os.environ.copy()

    install_spec = "-e.[test]"
    session.install(install_spec)
    session.run("pytest", *posargs, env=env)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def type_check(session):
    posargs = list(session.posargs)
    env = os.environ.copy()

    install_spec = ".[typing]"
    session.install(install_spec)
    session.run("mypy", "src", "tests", *posargs, env=env)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def ruff(session):
    posargs = list(session.posargs)
    env = os.environ.copy()

    install_spec = ".[lint]"
    session.install(install_spec)
    session.run("ruff", "check", *posargs, env=env)
