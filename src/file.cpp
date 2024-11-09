// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef _WIN32
// Workaround bug several symbol redefinition errors due to something including <winsock.h>
#include <winsock2.h>
#endif

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <hictk/balancing/methods.hpp>
#include <hictk/balancing/weights.hpp>
#include <hictk/bin_table.hpp>
#include <hictk/cooler/cooler.hpp>
#include <hictk/cooler/validation.hpp>
#include <hictk/file.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/hic.hpp>
#include <hictk/hic/common.hpp>
#include <hictk/hic/validation.hpp>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel_selector.hpp"
#include "hictkpy/reference.hpp"

namespace nb = nanobind;

namespace hictkpy::file {
static void ctor(hictk::File *fp, const std::filesystem::path &path,
                 std::optional<std::int32_t> resolution, std::string_view matrix_type,
                 std::string_view matrix_unit) {
  new (fp) hictk::File{path.string(), static_cast<std::uint32_t>(resolution.value_or(0)),
                       hictk::hic::ParseMatrixTypeStr(std::string{matrix_type}),
                       hictk::hic::ParseUnitStr(std::string{matrix_unit})};
}

static std::string repr(const hictk::File &f) {
  return fmt::format(FMT_STRING("File({})"), f.uri());
}

bool is_cooler(const std::filesystem::path &uri) {
  return bool(hictk::cooler::utils::is_cooler(uri.string()));
}

bool is_hic(const std::filesystem::path &uri) { return hictk::hic::utils::is_hic_file(uri); }

static hictkpy::PixelSelector fetch(const hictk::File &f, std::string_view range1,
                                    std::string_view range2, std::string_view normalization,
                                    std::string_view count_type, bool join,
                                    std::string_view query_type) {
  if (count_type != "float" && count_type != "int") {
    throw std::runtime_error(R"(count_type should be either "float" or "int")");
  }

  if (query_type != "UCSC" && query_type != "BED") {
    throw std::runtime_error("query_type should be either UCSC or BED");
  }

  if (normalization != "NONE") {
    count_type = "float";
  }

  if (range1.empty()) {
    assert(range2.empty());
    return std::visit(
        [&](const auto &ff) {
          auto sel = ff.fetch(hictk::balancing::Method{normalization});
          using SelT = decltype(sel);
          return hictkpy::PixelSelector(std::make_shared<const SelT>(std::move(sel)), count_type,
                                        join);
        },
        f.get());
  }

  if (range2.empty()) {
    range2 = range1;
  }

  const auto query_type_ =
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED;
  const auto gi1 = hictk::GenomicInterval::parse(f.chromosomes(), std::string{range1}, query_type_);
  const auto gi2 = hictk::GenomicInterval::parse(f.chromosomes(), std::string{range2}, query_type_);

  return std::visit(
      [&](const auto &ff) {
        // Workaround bug fixed in https://github.com/paulsengroup/hictk/pull/158
        auto sel = ff.fetch(fmt::format(FMT_STRING("{}"), gi1), fmt::format(FMT_STRING("{}"), gi2),
                            hictk::balancing::Method(normalization));

        using SelT = decltype(sel);
        return hictkpy::PixelSelector(std::make_shared<const SelT>(std::move(sel)), count_type,
                                      join);
      },
      f.get());
}

static nb::dict get_cooler_attrs(const hictk::cooler::File &clr) {
  nb::dict py_attrs;
  const auto &attrs = clr.attributes();

  py_attrs["bin-size"] = attrs.bin_size;
  py_attrs["bin-type"] = attrs.bin_type == hictk::BinTable::Type::fixed ? "fixed" : "variable";
  py_attrs["format"] = attrs.format;
  py_attrs["format-version"] = attrs.format_version;

  for (const auto &key : {"storage-mode", "creation-date", "generated-by", "assembly", "metadata",
                          "format-url", "nbins", "nchroms", "nnz", "sum", "cis"}) {
    py_attrs[key] = nb::none();
  }

  if (attrs.storage_mode.has_value()) {
    py_attrs["storage-mode"] = *attrs.storage_mode;
  }
  if (attrs.creation_date.has_value()) {
    py_attrs["creation-date"] = *attrs.creation_date;
  }
  if (attrs.generated_by.has_value()) {
    py_attrs["generated-by"] = *attrs.generated_by;
  }
  if (attrs.assembly.has_value()) {
    py_attrs["assembly"] = *attrs.assembly;
  }
  if (attrs.metadata.has_value()) {
    py_attrs["metadata"] = *attrs.metadata;
  }
  if (attrs.format_url.has_value()) {
    py_attrs["format-url"] = *attrs.format_url;
  }
  if (attrs.nbins.has_value()) {
    py_attrs["nbins"] = *attrs.nbins;
  }
  if (attrs.nchroms.has_value()) {
    py_attrs["nchroms"] = *attrs.nchroms;
  }
  if (attrs.nnz.has_value()) {
    py_attrs["nnz"] = *attrs.nnz;
  }
  if (attrs.sum.has_value()) {
    std::visit([&](const auto &sum) { py_attrs["sum"] = sum; }, *attrs.sum);
  }
  if (attrs.cis.has_value()) {
    std::visit([&](const auto &cis) { py_attrs["cis"] = cis; }, *attrs.cis);
  }

  return py_attrs;
}

static nb::dict get_hic_attrs(const hictk::hic::File &hf) {
  nb::dict py_attrs;

  py_attrs["bin_size"] = hf.resolution();
  py_attrs["format"] = "HIC";
  py_attrs["format_version"] = hf.version();
  py_attrs["assembly"] = hf.assembly();
  py_attrs["format-url"] = "https://github.com/aidenlab/hic-format";
  py_attrs["nbins"] = hf.bins().size();
  py_attrs["nchroms"] = hf.nchroms();

  return py_attrs;
}

static nb::dict attributes(const hictk::File &f) {
  if (f.is_cooler()) {
    return get_cooler_attrs(f.get<hictk::cooler::File>());
  }
  return get_hic_attrs(f.get<hictk::hic::File>());
}

static std::vector<std::string> avail_normalizations(const hictk::File &f) {
  const auto norms_ = f.avail_normalizations();
  std::vector<std::string> norms{norms_.size()};
  std::transform(norms_.begin(), norms_.end(), norms.begin(),
                 [](const auto &norm) { return norm.to_string(); });

  return norms;
}

static std::vector<double> weights(const hictk::File &f, std::string_view normalization,
                                   bool divisive) {
  const auto type = divisive ? hictk::balancing::Weights::Type::DIVISIVE
                             : hictk::balancing::Weights::Type::MULTIPLICATIVE;
  return f.normalization(normalization).to_vector(type);
}

static std::filesystem::path get_path(const hictk::File &f) { return f.path(); }

void declare_file_class(nb::module_ &m) {
  auto file = nb::class_<hictk::File>(m, "File",
                                      "Class representing a file handle to a .cool or .hic file.");

  file.def("__init__", &file::ctor, nb::arg("path"), nb::arg("resolution") = nb::none(),
           nb::arg("matrix_type") = "observed", nb::arg("matrix_unit") = "BP",
           "Construct a file object to a .hic, .cool or .mcool file given the file path and "
           "resolution.\n"
           "Resolution is ignored when opening single-resolution Cooler files.");

  file.def("__repr__", &file::repr, nb::rv_policy::take_ownership);

  file.def("uri", &hictk::File::uri, "Return the file URI.");
  file.def("path", &file::get_path, "Return the file path.", nb::rv_policy::take_ownership);

  file.def("is_hic", &hictk::File::is_hic, "Test whether file is in .hic format.");
  file.def("is_cooler", &hictk::File::is_cooler, "Test whether file is in .cool format.");

  file.def("chromosomes", &get_chromosomes_from_object<hictk::File>, nb::arg("include_ALL") = false,
           "Get chromosomes sizes as a dictionary mapping names to sizes.",
           nb::rv_policy::take_ownership);
  file.def("bins", &get_bins_from_object<hictk::File>, "Get table of bins.",
           nb::sig("def bins(self) -> hictkpy.BinTable"), nb::rv_policy::take_ownership);

  file.def("resolution", &hictk::File::resolution, "Get the bin size in bp.");
  file.def("nbins", &hictk::File::nbins, "Get the total number of bins.");
  file.def("nchroms", &hictk::File::nchroms, "Get the total number of chromosomes.");

  file.def("attributes", &file::attributes, "Get file attributes as a dictionary.",
           nb::rv_policy::take_ownership);

  file.def("fetch", &file::fetch, nb::keep_alive<0, 1>(), nb::arg("range1") = "",
           nb::arg("range2") = "", nb::arg("normalization") = "NONE", nb::arg("count_type") = "int",
           nb::arg("join") = false, nb::arg("query_type") = "UCSC",
           "Fetch interactions overlapping a region of interest.", nb::rv_policy::take_ownership);

  file.def("avail_normalizations", &file::avail_normalizations,
           "Get the list of available normalizations.", nb::rv_policy::take_ownership);
  file.def("has_normalization", &hictk::File::has_normalization, nb::arg("normalization"),
           "Check whether a given normalization is available.");
  file.def("weights", &file::weights, nb::arg("name"), nb::arg("divisive") = true,
           "Fetch the balancing weights for the given normalization method.");
}

}  // namespace hictkpy::file
