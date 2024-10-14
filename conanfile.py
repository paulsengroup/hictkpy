# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.tools.build import check_min_cppstd

required_conan_version = ">=1.53.0"


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

    @property
    def _min_cppstd(self):
        return 17

    def requirements(self):
        self.requires("bshoshany-thread-pool/4.1.0#be1802a8768416a6c9b1393cf0ce5e9c")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("eigen/3.4.0#2e192482a8acff96fe34766adca2b24c")
        self.requires("fast_float/6.1.5#e067b96a6271d1b4c255858ca9805bdd")
        self.requires("fmt/11.0.2#5c7438ef4d5d69ab106a41e460ce11f3", force=True)
        self.requires("hdf5/1.14.4.3#df1467d7374938c231edbe10e83f2bb4", force=True)
        self.requires("highfive/2.10.0#3d1bd25944a57fa1bc30a0a22923d528")
        self.requires("libdeflate/1.22#f95aebe763153ccbc4cc76c023e42e5a")
        self.requires("parallel-hashmap/1.4.0#36ac84df77219748440cdb0f23624d56")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")
        self.requires("spdlog/1.14.1#972bbf70be1da4bc57ea589af0efde03")
        self.requires("zstd/1.5.6#afefe79a309bc2a7b9f56c2093504c8b", force=True)

    def validate(self):
        if self.settings.get_safe("compiler.cppstd"):
            check_min_cppstd(self, self._min_cppstd)

    def configure(self):
        if self.settings.compiler in ["clang", "gcc"]:
            self.settings.compiler.libcxx = "libstdc++11"

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
