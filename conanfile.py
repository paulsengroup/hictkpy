# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.tools.microsoft import is_msvc

required_conan_version = ">=2"


class HictkpyConan(ConanFile):
    name = "hictkpy"
    description = "Python bindings for hictk: read and write .cool and .hic files directly from Python."
    license = "MIT"
    topics = ("hictk", "hictkpy", "bioinformatics")
    homepage = "https://github.com/paulsengroup/hictkpy"
    url = "https://github.com/paulsengroup/hictkpy"
    package_type = "shared-library"
    settings = "os", "arch", "compiler", "build_type"

    options = {
        "shared": [True, False],
        "fPIC": [True, False],
    }

    default_options = {
        "shared": False,
        "fPIC": True,
    }

    generators = "CMakeDeps"

    def requirements(self):
        # Arrow 21.0.0 can't find certain kernels (e.g., sort_indices)
        self.requires("arrow/20.0.0#6e04404a336dd16f08062f6923e6f8f1")
        if is_msvc(self):
            self.requires("boost/1.88.0#14ecfc01dd5a690f15e1318e56a6b78c", force=True)  # arrow
        self.requires("bshoshany-thread-pool/5.0.0#d94da300363f0c35b8f41b2c5490c94d")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("eigen/3.4.90-unstable+git.2025.08.15#b407f03f085cdb246f6bcbadd84fe9db", force=True)
        self.requires("fast_float/8.0.2#846ad0ebab16bc265c511095c3b490e9")
        self.requires("fmt/11.2.0#579bb2cdf4a7607621beea4eb4651e0f", force=True)
        self.requires("hdf5/1.14.6#a5cdabab5e051941b7c7a220005d1d60", force=True)
        self.requires("highfive/2.10.0#75c849a0d940b2d4dae6055915132690")
        self.requires("libdeflate/1.23#4994bea7cf7e93789da161fac8e26a53")
        self.requires("parallel-hashmap/2.0.0#82acae64ffe2693fff5fb3f9df8e1746")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")
        self.requires("spdlog/1.15.3#3ca0e9e6b83af4d0151e26541d140c86")
        self.requires("zstd/1.5.7#fde461c0d847a22f16d3066774f61b11", force=True)

    def configure(self):
        self.options["arrow"].compute = True
        self.options["arrow"].filesystem_layer = False
        self.options["arrow"].parquet = False
        self.options["arrow"].with_boost = is_msvc(self)
        self.options["arrow"].with_re2 = True
        self.options["arrow"].with_thrift = False
        if is_msvc(self):
            self.options["boost"].system_no_deprecated = True
            self.options["boost"].asio_no_deprecated = True
            self.options["boost"].filesystem_no_deprecated = True
            self.options["boost"].filesystem_version = 4
            self.options["boost"].zlib = False
            self.options["boost"].bzip2 = False
            self.options["boost"].lzma = False
            self.options["boost"].zstd = False
            self.options["boost"].without_atomic = False
            self.options["boost"].without_charconv = True
            self.options["boost"].without_chrono = True
            self.options["boost"].without_cobalt = True
            self.options["boost"].without_container = True
            self.options["boost"].without_context = True
            self.options["boost"].without_contract = True
            self.options["boost"].without_coroutine = True
            self.options["boost"].without_date_time = True
            self.options["boost"].without_exception = True
            self.options["boost"].without_fiber = True
            self.options["boost"].without_filesystem = False
            self.options["boost"].without_graph = True
            self.options["boost"].without_graph_parallel = True
            self.options["boost"].without_iostreams = True
            self.options["boost"].without_json = True
            self.options["boost"].without_locale = True
            self.options["boost"].without_log = True
            self.options["boost"].without_math = True
            self.options["boost"].without_mpi = True
            self.options["boost"].without_nowide = True
            self.options["boost"].without_process = True
            self.options["boost"].without_program_options = True
            self.options["boost"].without_python = True
            self.options["boost"].without_random = True
            self.options["boost"].without_regex = True
            self.options["boost"].without_serialization = True
            self.options["boost"].without_stacktrace = True
            self.options["boost"].without_system = False
            self.options["boost"].without_test = True
            self.options["boost"].without_thread = True
            self.options["boost"].without_timer = True
            self.options["boost"].without_type_erasure = True
            self.options["boost"].without_url = True
            self.options["boost"].without_wave = True
        self.options["fmt"].header_only = True
        self.options["hdf5"].enable_cxx = False
        self.options["hdf5"].hl = False
        self.options["hdf5"].threadsafe = False
        self.options["hdf5"].parallel = False
        self.options["highfive"].with_boost = False
        self.options["highfive"].with_eigen = False
        self.options["highfive"].with_opencv = False
        self.options["highfive"].with_xtensor = False
        self.options["spdlog"].header_only = True
        self.options["zstd"].build_programs = False
