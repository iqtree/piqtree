"""Decorators for IQ-TREE functions."""

import io
import os
import pathlib
import sys
import tempfile
from collections.abc import Callable
from functools import wraps

from piqtree.exceptions import IqTreeError


def _fd_or_fallback(stream: object, fallback_fd: int) -> int:
    """Return the file descriptor for the stream, or a fallback file descriptor.

    Parameters
    ----------
    stream : stdout or stderr stream
        The io stream.
    fallback_fd : int
        Fallback file descriptor.

    Returns
    -------
    int
        The file descriptor for the stream, or the fallback on error.
    """
    try:
        return stream.fileno()  # type: ignore[attr-defined]
    except io.UnsupportedOperation:  # pragma: no cover
        return fallback_fd  # pragma: no cover


def iqtree_func[**Param, RetType](
    func: Callable[Param, RetType],
    *,
    hide_files: bool = False,
) -> Callable[Param, RetType]:
    """IQ-TREE function wrapper.

    Hides stdout and stderr, as well as any output files.

    Parameters
    ----------
    func : Callable[Param, RetType]
        The IQ-TREE library function.
    hide_files : bool, optional
        Whether hiding output files is necessary, by default False.

    Returns
    -------
    Callable[Param, RetType]
        The wrapped IQ-TREE function.

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

        # Fetch file descriptors
        out_fd = _fd_or_fallback(sys.stdout, 1)
        err_fd = _fd_or_fallback(sys.stderr, 2)

        # Save original stdout and stderr file descriptors
        saved_stdout_fd = os.dup(out_fd)
        saved_stderr_fd = os.dup(err_fd)

        # Open /dev/null (or NUL on Windows) as destination for stdout and stderr
        devnull_fd = os.open(os.devnull, os.O_WRONLY)

        if hide_files:
            original_dir = pathlib.Path.cwd()
            tempdir = tempfile.TemporaryDirectory(prefix=f"piqtree_{func.__name__}")
            os.chdir(tempdir.name)

        try:
            # Replace stdout and stderr with /dev/null
            os.dup2(devnull_fd, out_fd)
            os.dup2(devnull_fd, err_fd)

            # Call the wrapped function
            return func(*args, **kwargs)
        except RuntimeError as e:
            raise IqTreeError(e) from None
        finally:
            # Flush stdout and stderr
            sys.stdout.flush()
            sys.stderr.flush()

            # Restore stdout and stderr
            os.dup2(saved_stdout_fd, out_fd)
            os.dup2(saved_stderr_fd, err_fd)

            # Close the devnull file descriptor, and duplicated file descriptors
            os.close(devnull_fd)
            os.close(saved_stdout_fd)
            os.close(saved_stderr_fd)

            if hide_files:
                os.chdir(original_dir)
                tempdir.cleanup()

    return wrapper_iqtree_func
