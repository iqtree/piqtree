"""setup for piqtree."""

import os
import platform
import subprocess
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

LIBRARY_DIR = "src/piqtree/_libiqtree/"
IQTREE_LIB_NAME = "iqtree2"


def get_brew_prefix(package: str) -> Path:
    """Get the prefix path for a specific Homebrew package."""
    return Path(
        subprocess.check_output(["brew", "--prefix", package]).strip().decode("utf-8"),
    )


extra_libs = []
extra_compile_args = []
include_dirs = []
library_dirs = []


def include_dlls() -> None:
    import shutil

    from delvewheel._dll_utils import get_all_needed

    needed_dll_paths, _, _, _ = get_all_needed(
        LIBRARY_DIR + f"{IQTREE_LIB_NAME}.dll",
        set(),
        None,
        "raise",
        False,  # noqa: FBT003
        False,  # noqa: FBT003
        0,
    )

    for dll_path in needed_dll_paths:
        shutil.copy(dll_path, LIBRARY_DIR)


def setup_windows() -> None:
    include_dlls()


def setup_macos() -> None:
    brew_prefix_llvm = get_brew_prefix("llvm")
    brew_prefix_libomp = get_brew_prefix("libomp")

    # Use Homebrew's clang/clang++
    os.environ["CC"] = str(brew_prefix_llvm / "bin" / "clang")
    os.environ["CXX"] = str(brew_prefix_llvm / "bin" / "clang++")

    # Define OpenMP flags and libraries for macOS
    extra_compile_args.extend(["-Xpreprocessor", "-fopenmp"])
    extra_libs.extend(["z", "omp"])

    # Use the paths from Homebrew for libomp
    include_dirs.extend([str(brew_prefix_libomp / "include")])
    library_dirs.extend(
        [
            str(brew_prefix_libomp / "lib"),
            str(brew_prefix_llvm / "lib"),
        ],
    )


def setup_linux() -> None:
    extra_compile_args.extend(["-fopenmp"])
    extra_libs.extend(["z", "gomp"])


match system := platform.system():
    case "Windows":
        setup_windows()
    case "Darwin":
        setup_macos()
    case "Linux":
        setup_linux()
    case _:
        msg = f"Unsupported platform: {system}"
        raise ValueError(msg)

ext_modules = [
    Pybind11Extension(
        "_piqtree",
        ["src/piqtree/_libiqtree/_piqtree.cpp"],
        library_dirs=[
            *library_dirs,
            LIBRARY_DIR,
        ],
        libraries=[IQTREE_LIB_NAME, *extra_libs],
        extra_compile_args=extra_compile_args,
        include_dirs=include_dirs,
    ),
]

setup(
    name="piqtree",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    package_data={"piqtree": ["_libiqtree/*.dll"]},
)
