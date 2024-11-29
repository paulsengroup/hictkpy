// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/multires_file.hpp"

#include <fmt/format.h>

#include <hictk/cooler/validation.hpp>
#include <hictk/multires_file.hpp>
#include <string>
#include <vector>

#include "hictkpy/nanobind.hpp"
#include "hictkpy/reference.hpp"

namespace nb = nanobind;

namespace hictkpy::multires_file {

static void ctor(hictk::MultiResFile* fp, const std::filesystem::path& path) {
  new (fp) hictk::MultiResFile{path.string()};
}

static std::string repr(const hictk::MultiResFile& mrf) {
  return fmt::format(FMT_STRING("MultiResFile({})"), mrf.path());
}

static std::filesystem::path get_path(const hictk::MultiResFile& mrf) { return mrf.path(); }

[[nodiscard]] static auto get_resolutions(const hictk::MultiResFile& f) {
  using WeightVector = nb::ndarray<nb::numpy, nb::shape<-1>, nb::c_contig, std::uint32_t>;

  // NOLINTNEXTLINE
  auto* resolutions_ptr = new std::vector<std::uint32_t>(f.resolutions());

  auto capsule = nb::capsule(resolutions_ptr, [](void* vect_ptr) noexcept {
    delete reinterpret_cast<std::vector<std::uint32_t>*>(vect_ptr);  // NOLINT
  });

  return WeightVector{resolutions_ptr->data(), {resolutions_ptr->size()}, capsule};
}

static nb::dict get_attrs(const hictk::hic::File& hf) {
  nb::dict py_attrs;

  py_attrs["format"] = "HIC";
  py_attrs["format-version"] = hf.version();
  py_attrs["assembly"] = hf.assembly();
  py_attrs["format-url"] = "https://github.com/aidenlab/hic-format";
  py_attrs["nchroms"] = hf.nchroms();

  for (const auto& [k, v] : hf.attributes()) {
    py_attrs[nb::cast(k)] = v;
  }

  return py_attrs;
}

static nb::dict get_attrs(const hictk::cooler::MultiResFile& mclr) {
  nb::dict py_attrs;

  py_attrs["format"] = mclr.attributes().format;
  py_attrs["format-version"] = mclr.attributes().format_version;
  py_attrs["format-url"] = "https://github.com/open2c/cooler";
  py_attrs["assembly"] =
      mclr.open(mclr.resolutions().front()).attributes().assembly.value_or("unknown");
  py_attrs["nchroms"] = mclr.chromosomes().size();

  return py_attrs;
}

static nb::dict attributes(const hictk::MultiResFile& f) {
  auto attrs = f.is_hic() ? get_attrs(f.open(f.resolutions().front()).get<hictk::hic::File>())
                          : get_attrs(hictk::cooler::MultiResFile{f.path()});
  attrs["resolutions"] = get_resolutions(f);

  return attrs;
}

bool is_mcool_file(const std::filesystem::path& path) {
  return bool(hictk::cooler::utils::is_multires_file(path.string()));
}

void declare_multires_file_class(nb::module_& m) {
  auto mres_file = nb::class_<hictk::MultiResFile>(
      m, "MultiResFile", "Class representing a file handle to a .hic or .mcool file");
  mres_file.def("__init__", &multires_file::ctor, nb::arg("path"),
                "Open a multi-resolution Cooler file (.mcool) or .hic file.");

  mres_file.def("__repr__", &multires_file::repr, nb::rv_policy::move);

  mres_file.def("path", &multires_file::get_path, "Get the file path.", nb::rv_policy::move);
  mres_file.def("is_mcool", &hictk::MultiResFile::is_mcool,
                "Test whether the file is in .mcool format.");
  mres_file.def("is_hic", &hictk::MultiResFile::is_hic, "Test whether the file is in .hic format.");
  mres_file.def("chromosomes", &get_chromosomes_from_object<hictk::MultiResFile>,
                nb::arg("include_ALL") = false,
                "Get chromosomes sizes as a dictionary mapping names to sizes.",
                nb::rv_policy::take_ownership);
  mres_file.def("resolutions", &get_resolutions, "Get the list of available resolutions.",
                nb::rv_policy::take_ownership);
  mres_file.def("attributes", &multires_file::attributes, "Get file attributes as a dictionary.",
                nb::rv_policy::take_ownership);
  mres_file.def("__getitem__", &hictk::MultiResFile::open,
                "Open the Cooler or .hic file corresponding to the resolution given as input.",
                nb::rv_policy::move);
}

}  // namespace hictkpy::multires_file
