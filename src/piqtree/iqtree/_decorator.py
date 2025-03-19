"""Decorators for IQ-TREE functions."""

import os
import pathlib
import sys
import tempfile
from collections.abc import Callable
from functools import wraps
from typing import TypeVar

from typing_extensions import ParamSpec

from piqtree.exceptions import IqTreeError

Param = ParamSpec("Param")
RetType = TypeVar("RetType")


def iqtree_func(
    func: Callable[Param, RetType],
    *,
    hide_files: bool | None = False,
) -> Callable[Param, RetType]:
    """IQ-TREE function wrapper.

    Hides stdout and stderr, as well as any output files.

    Parameters
    ----------
    func : Callable[Param, RetType]
        The IQ-TREE library function.
    hide_files : bool | None, optional
        Whether hiding output files is necessary, by default False.

    Returns
    -------
    Callable[Param, RetType]
        The wrappe IQ-TREE function.

    Raises
    ------
    IqTreeError
        An error from the IQ-TREE library.

    """

    @wraps(func)
    def wrapper_iqtree_func(*args: Param.args, **kwargs: Param.kwargs) -> RetType:
        # Flush stdout and stderr
        sys.stdout.flush()
        sys.stderr.flush()

        # Save original stdout and stderr file descriptors
        original_stdout_fd = os.dup(sys.stdout.fileno())
        original_stderr_fd = os.dup(sys.stderr.fileno())

        # Open /dev/null (or NUL on Windows) as destination for stdout and stderr
        devnull_fd = os.open(os.devnull, os.O_WRONLY)

        if hide_files:
            original_dir = pathlib.Path.cwd()
            tempdir = tempfile.TemporaryDirectory(prefix=f"piqtree_{func.__name__}")
            os.chdir(tempdir.name)

        try:
            # Replace stdout and stderr with /dev/null
            os.dup2(devnull_fd, sys.stdout.fileno())
            os.dup2(devnull_fd, sys.stderr.fileno())

            # Call the wrapped function
            return func(*args, **kwargs)
        except RuntimeError as e:
            raise IqTreeError(e) from None
        finally:
            # Flush stdout and stderr
            sys.stdout.flush()
            sys.stderr.flush()

            # Restore stdout and stderr
            os.dup2(original_stdout_fd, sys.stdout.fileno())
            os.dup2(original_stderr_fd, sys.stderr.fileno())

            # Close the devnull file descriptor
            os.close(devnull_fd)

            if hide_files:
                os.chdir(original_dir)
                tempdir.cleanup()

    return wrapper_iqtree_func
