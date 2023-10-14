#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext: CMakeExtension) -> None:
        # Must be in this form due to bug in .resolve() only fixed in Python 3.10+
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        extdir = ext_fullpath.parent.resolve() / ext.name
        extdir.mkdir(exist_ok=True, parents=True)
        with open(extdir / "__init__.py", "w") as f:
            symbols = ", ".join(("File", "PixelSelector", "is_cooler", "is_hic", "__hictk_version__"))
            f.write(
                textwrap.dedent(
                    f"""
                    from .hictkpy import {symbols}
                    from .hictkpy import cooler

                    try:
                        from importlib.metadata import version
                    except ModuleNotFoundError:
                        from importlib_metadata import version

                    __version__ = version("hictkpy")
                    """
                )
            )

        # Using this requires trailing slash for auto-detection & inclusion of
        # auxiliary "native" libs

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        cmake_args = [
            "-DHICTK_ENABLE_TESTING=OFF",
            "-DHICTK_BUILD_EXAMPLES=OFF",
            "-DHICTK_BUILD_BENCHMARKS=OFF",
            "-DHICTK_BUILD_TOOLS=OFF",
            "-DBUILD_SHARED_LIBS=ON",
            "-DHICTK_ENABLE_GIT_VERSION_TRACKING=OFF",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
        ]
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja

                    ninja_executable_path = Path(ninja.BIN_DIR) / "ninja"
                    cmake_args += [
                        "-GNinja",
                        f"-DCMAKE_MAKE_PROGRAM:FILEPATH={ninja_executable_path}",
                    ]
                except ImportError:
                    pass

        else:
            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        if "HICTKPY_SETUP_SKIP_CONAN" not in os.environ:
            # profile detect fails if profile already exists
            subprocess.run(["conan", "profile", "detect"], check=False)
            subprocess.run(["conan", "profile", "detect", "--name", self.plat_name], check=False)
            subprocess.run(
                [
                    "conan",
                    "install",
                    ext.sourcedir,
                    f"-pr:b=default",
                    f"-pr:h={self.plat_name}",
                    "-s",
                    f"build_type={cfg}",
                    "-s",
                    "compiler.cppstd=17",
                    f"--output-folder={build_temp.absolute()}",
                    "-o",
                    "*/*:shared=True",
                    "--build=missing",
                ],
                check=True,
            )
            cmake_args += [f"-DCMAKE_PREFIX_PATH={build_temp.absolute()}"]

        subprocess.run(["cmake", ext.sourcedir, *cmake_args], cwd=build_temp, check=True)
        subprocess.run(["cmake", "--build", ".", *build_args], cwd=build_temp, check=True)


# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    ext_modules=[CMakeExtension("hictkpy")],
    cmdclass={"build_ext": CMakeBuild},
)
