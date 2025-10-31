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

    def build_requirements(self):
        # re2
        self.requires("abseil/20250814.0#4e0fdd34a888b97aca482e648fc27a3b", force=True)
        if is_msvc(self):
            # arrow
            self.requires("boost/1.89.0#010f59feedfd171b15f467b42f723d13", force=True)
        # arrow
        self.requires("re2/20250722#7547baba4648ebb432652af97ec9c972", force=True)
        # hdf5
        self.requires("zlib/1.3.1#b8bc2603263cf7eccbd6e17e66b0ed76", force=True)

    def requirements(self):
        self.requires("arrow/22.0.0#0a31a3dca837570a86d22b65ded024f4")
        self.requires("bshoshany-thread-pool/5.0.0#d94da300363f0c35b8f41b2c5490c94d")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("eigen/5.0.0#f7561f543f4aafd6d2dc1f6d677e3075", force=True)
        self.requires("fast_float/8.1.0#bbf67486bec084d167da0f3e13eee534")
        self.requires("fmt/12.1.0#50abab23274d56bb8f42c94b3b9a40c7", force=True)
        self.requires("hdf5/1.14.6#0b780319690d537e6cb0683244919955", force=True)
        self.requires("highfive/3.1.1#d0c724526ebc8ce396ffa1bf7f3c7b64")
        self.requires("libdeflate/1.23#4994bea7cf7e93789da161fac8e26a53")
        self.requires("parallel-hashmap/2.0.0#82acae64ffe2693fff5fb3f9df8e1746")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")
        self.requires("spdlog/1.16.0#942c2c39562ae25ba575d9c8e2bdf3b6")
        self.requires("zstd/1.5.7#fde461c0d847a22f16d3066774f61b11", force=True)

    def configure(self):
        self.options["arrow"].compute = True
        self.options["arrow"].filesystem_layer = False
        self.options["arrow"].parquet = False
        self.options["arrow"].with_boost = is_msvc(self)
        self.options["arrow"].with_re2 = True
        self.options["arrow"].with_thrift = False
        if is_msvc(self):
            self.options["boost"].header_only = True
            self.options["boost"].system_no_deprecated = True
            self.options["boost"].asio_no_deprecated = True
            self.options["boost"].filesystem_no_deprecated = True
            self.options["boost"].filesystem_version = 4
        self.options["eigen"].MPL2_only = True
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
