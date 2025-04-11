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
#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel_selector.hpp"
#include "hictkpy/reference.hpp"
#include "hictkpy/to_pyarrow.hpp"

namespace nb = nanobind;

namespace hictkpy::file {
static void ctor(hictk::File *fp, const std::filesystem::path &path,
                 std::optional<std::int32_t> resolution, std::string_view matrix_type,
                 std::string_view matrix_unit) {
  if (resolution.value_or(0) < 0) {
    throw std::invalid_argument("resolution must be non-negative");
  }

  std::optional<std::uint32_t> resolution_{};
  if (resolution.has_value()) {
    resolution_ = static_cast<std::uint32_t>(*resolution);
  }

  new (fp) hictk::File{path.string(), resolution_,
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

static hictkpy::PixelSelector fetch(const hictk::File &f, std::optional<std::string_view> range1,
                                    std::optional<std::string_view> range2,
                                    std::optional<std::string_view> normalization,
                                    std::string_view count_type, bool join,
                                    std::string_view query_type,
                                    std::optional<std::int64_t> diagonal_band_width) {
  if (count_type != "float" && count_type != "int") {
    throw std::runtime_error(R"(count_type should be either "float" or "int")");
  }

  if (query_type != "UCSC" && query_type != "BED") {
    throw std::runtime_error("query_type should be either UCSC or BED");
  }

  const hictk::balancing::Method normalization_method{normalization.value_or("NONE")};

  if (normalization_method != hictk::balancing::Method::NONE()) {
    count_type = "float";
  }

  if (!range1.has_value() || range1->empty()) {
    assert(!range2.has_value() || range2->empty());
    return std::visit(
        [&](const auto &ff) {
          using FileT = remove_cvref_t<decltype(ff)>;
          if constexpr (std::is_same_v<FileT, hictk::cooler::File>) {
            auto sel = ff.fetch(normalization_method, diagonal_band_width.has_value());
            using SelT = decltype(sel);
            return hictkpy::PixelSelector(std::make_shared<const SelT>(std::move(sel)), count_type,
                                          join, diagonal_band_width);
          } else {
            auto sel = ff.fetch(normalization_method, diagonal_band_width);
            using SelT = decltype(sel);
            return hictkpy::PixelSelector(std::make_shared<const SelT>(std::move(sel)), count_type,
                                          join, diagonal_band_width);
          }
        },
        f.get());
  }

  if (!range2.has_value() || range2->empty()) {
    range2 = range1;
  }

  const auto query_type_ =
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED;
  const auto gi1 =
      hictk::GenomicInterval::parse(f.chromosomes(), std::string{*range1}, query_type_);
  const auto gi2 =
      hictk::GenomicInterval::parse(f.chromosomes(), std::string{*range2}, query_type_);

  return std::visit(
      [&](const auto &ff) {
        using FileT = remove_cvref_t<decltype(ff)>;
        if constexpr (std::is_same_v<FileT, hictk::hic::File>) {
          auto sel = ff.fetch(gi1.chrom().name(), gi1.start(), gi1.end(), gi2.chrom().name(),
                              gi2.start(), gi2.end(), normalization_method, diagonal_band_width);

          using SelT = decltype(sel);
          return hictkpy::PixelSelector(std::make_shared<const SelT>(std::move(sel)), count_type,
                                        join, diagonal_band_width);
        } else {
          auto sel = ff.fetch(gi1.chrom().name(), gi1.start(), gi1.end(), gi2.chrom().name(),
                              gi2.start(), gi2.end(), normalization_method);

          using SelT = decltype(sel);
          return hictkpy::PixelSelector(std::make_shared<const SelT>(std::move(sel)), count_type,
                                        join, diagonal_band_width);
        }
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

  py_attrs["bin-size"] = hf.resolution();
  py_attrs["format"] = "HIC";
  py_attrs["format-version"] = hf.version();
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

static auto weights(const hictk::File &f, std::string_view normalization, bool divisive) {
  using WeightVector = nb::ndarray<nb::numpy, nb::shape<-1>, nb::c_contig, double>;

  if (normalization == "NONE") {
    return WeightVector{};
  }

  const auto type = divisive ? hictk::balancing::Weights::Type::DIVISIVE
                             : hictk::balancing::Weights::Type::MULTIPLICATIVE;

  // NOLINTNEXTLINE
  auto *weights_ptr = new std::vector<double>(f.normalization(normalization).to_vector(type));

  auto capsule = nb::capsule(weights_ptr, [](void *vect_ptr) noexcept {
    delete reinterpret_cast<std::vector<double> *>(vect_ptr);  // NOLINT
  });

  return WeightVector{weights_ptr->data(), {weights_ptr->size()}, capsule};
}

static nb::object weights_df(const hictk::File &f, const std::vector<std::string> &normalizations,
                             bool divisive) {
  phmap::flat_hash_set<std::string_view> names(normalizations.size());
  arrow::FieldVector fields(normalizations.size());
  std::vector<std::shared_ptr<arrow::Array>> columns(normalizations.size());

  fields.clear();
  columns.clear();

  const auto type = divisive ? hictk::balancing::Weights::Type::DIVISIVE
                             : hictk::balancing::Weights::Type::MULTIPLICATIVE;

  for (const auto &normalization : normalizations) {
    if (normalization == "NONE") {
      fields.resize(fields.size() - 1);
      columns.resize(columns.size() - 1);
      continue;
    }

    if (names.contains(normalization)) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("found duplicated value \"{}\" in the provided normalization name list"),
          normalization));
    }

    names.emplace(normalization);
    fields.emplace_back(arrow::field(normalization, arrow::float64(), false));
    columns.emplace_back(std::make_shared<arrow::DoubleArray>(
        f.nbins(), arrow::Buffer::FromVector(f.normalization(normalization).to_vector(type)),
        nullptr, 0, 0));
  }

  auto schema = std::make_shared<arrow::Schema>(std::move(fields));

  return export_pyarrow_table(arrow::Table::Make(std::move(schema), columns))
      .attr("to_pandas")(nb::arg("self_destruct") = true);
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

  file.def("__repr__", &file::repr, nb::rv_policy::move);

  file.def("uri", &hictk::File::uri, "Return the file URI.", nb::rv_policy::move);
  file.def("path", &file::get_path, "Return the file path.", nb::rv_policy::move);

  file.def("is_hic", &hictk::File::is_hic, "Test whether file is in .hic format.");
  file.def("is_cooler", &hictk::File::is_cooler, "Test whether file is in .cool format.");

  file.def("chromosomes", &get_chromosomes_from_object<hictk::File>, nb::arg("include_ALL") = false,
           "Get chromosomes sizes as a dictionary mapping names to sizes.",
           nb::rv_policy::take_ownership);
  file.def("bins", &get_bins_from_object<hictk::File>, "Get table of bins.",
           nb::sig("def bins(self) -> hictkpy.BinTable"), nb::rv_policy::move);

  file.def("resolution", &hictk::File::resolution, "Get the bin size in bp.");
  file.def("nbins", &hictk::File::nbins, "Get the total number of bins.");
  file.def("nchroms", &hictk::File::nchroms, nb::arg("include_ALL") = false,
           "Get the total number of chromosomes.");

  file.def("attributes", &file::attributes, "Get file attributes as a dictionary.",
           nb::rv_policy::take_ownership);

  file.def("fetch", &file::fetch, nb::keep_alive<0, 1>(), nb::arg("range1") = nb::none(),
           nb::arg("range2") = nb::none(), nb::arg("normalization") = nb::none(),
           nb::arg("count_type") = "int", nb::arg("join") = false, nb::arg("query_type") = "UCSC",
           nb::arg("diagonal_band_width") = nb::none(),
           "Fetch interactions overlapping a region of interest.", nb::rv_policy::move);

  file.def("avail_normalizations", &file::avail_normalizations,
           "Get the list of available normalizations.", nb::rv_policy::move);
  file.def("has_normalization", &hictk::File::has_normalization, nb::arg("normalization"),
           "Check whether a given normalization is available.");
  file.def("weights", &file::weights, nb::arg("name"), nb::arg("divisive") = true,
           "Fetch the balancing weights for the given normalization method.",
           nb::sig("def weights(self, name: str, divisive: bool = True) -> "
                   "numpy.ndarray[float]"),
           nb::rv_policy::take_ownership);
  file.def(
      "weights", &file::weights_df, nb::arg("names"), nb::arg("divisive") = true,
      "Fetch the balancing weights for the given normalization methods."
      "Weights are returned as a pandas.DataFrame.",
      nb::sig("def weights(self, names: collections.abc.Sequence[str], divisive: bool = True) -> "
              "pandas.DataFrame"),
      nb::rv_policy::take_ownership);
}

}  // namespace hictkpy::file
