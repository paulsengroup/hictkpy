// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/file.hpp"

#ifdef _WIN32
// Workaround bug several symbol redefinition errors due to something including <winsock.h>
#include <winsock2.h>
#endif

#include <H5public.h>
#include <arrow/array/array_base.h>
#include <arrow/buffer.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

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
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel_selector.hpp"
#include "hictkpy/reference.hpp"
#include "hictkpy/table.hpp"
#include "hictkpy/to_numpy.hpp"
#include "hictkpy/type.hpp"
#include "hictkpy/variant.hpp"

namespace nb = nanobind;

namespace hictkpy {
static hictkpy::PixelSelector fetch_gw_impl(const hictk::File &f,
                                            const hictk::balancing::Method &normalization,
                                            NumericDtype count_type, bool join,
                                            std::optional<std::int64_t> diagonal_band_width) {
  return std::visit(
      [&](const auto &ff) -> hictkpy::PixelSelector {
        using FileT = remove_cvref_t<decltype(ff)>;
        if constexpr (std::is_same_v<FileT, hictk::cooler::File>) {
          auto sel = [&]() {
            HICTKPY_LOCK_COOLER_MTX_SCOPED
            return ff.fetch(normalization, diagonal_band_width.has_value());
          }();
          using SelT = decltype(sel);
          return {std::make_shared<const SelT>(std::move(sel)), count_type, join,
                  diagonal_band_width};
        } else {
          auto sel = ff.fetch(normalization, diagonal_band_width);
          using SelT = decltype(sel);
          return {std::make_shared<const SelT>(std::move(sel)), count_type, join,
                  diagonal_band_width};
        }
      },
      f.get());
}

static hictkpy::PixelSelector fetch_impl(const hictk::File &f, const hictk::GenomicInterval &range1,
                                         const hictk::GenomicInterval &range2,
                                         const hictk::balancing::Method &normalization,
                                         NumericDtype count_type, bool join,
                                         std::optional<std::int64_t> diagonal_band_width) {
  const auto chrom1 = range1.chrom().name();
  const auto start1 = range1.start();
  const auto end1 = range1.end();
  const auto chrom2 = range2.chrom().name();
  const auto start2 = range2.start();
  const auto end2 = range2.end();

  return std::visit(
      [&](const auto &ff) -> hictkpy::PixelSelector {
        using FileT = remove_cvref_t<decltype(ff)>;
        if constexpr (std::is_same_v<FileT, hictk::hic::File>) {
          auto sel = ff.fetch(chrom1, start1, end1, chrom2, start2, end2, normalization,
                              diagonal_band_width);

          using SelT = decltype(sel);
          return {std::make_shared<const SelT>(std::move(sel)), count_type, join,
                  diagonal_band_width};
        } else {
          auto sel = [&]() {
            HICTKPY_LOCK_COOLER_MTX_SCOPED
            return ff.fetch(chrom1, start1, end1, chrom2, start2, end2, normalization);
          }();

          using SelT = decltype(sel);
          return {std::make_shared<const SelT>(std::move(sel)), count_type, join,
                  diagonal_band_width};
        }
      },
      f.get());
}

static hictkpy::PixelSelector fetch_impl(const hictk::File &f,
                                         const std::optional<std::string_view> &range1,
                                         const std::optional<std::string_view> &range2,
                                         const hictk::balancing::Method &normalization,
                                         NumericDtype count_type, bool join,
                                         hictk::GenomicInterval::Type query_type,
                                         std::optional<std::int64_t> diagonal_band_width) {
  if (!range1.has_value() || range1->empty()) {
    assert(!range2.has_value() || range2->empty());
    return fetch_gw_impl(f, normalization, count_type, join, diagonal_band_width);
  }

  if (!range2.has_value() || range2->empty()) {
    return fetch_impl(f, range1, range1, normalization, count_type, join, query_type,
                      diagonal_band_width);
  }

  return fetch_impl(
      f, hictk::GenomicInterval::parse(f.chromosomes(), std::string{*range1}, query_type),
      hictk::GenomicInterval::parse(f.chromosomes(), std::string{*range2}, query_type),
      normalization, count_type, join, diagonal_band_width

  );
}

static hictkpy::PixelSelector fetch(const File &f, const std::optional<std::string_view> &range1,
                                    const std::optional<std::string_view> &range2,
                                    const std::optional<std::string_view> &normalization,
                                    std::variant<nb::type_object, std::string_view> count_type,
                                    bool join, std::string_view query_type,
                                    std::optional<std::int64_t> diagonal_band_width) {
  if (query_type != "UCSC" && query_type != "BED") {
    throw std::runtime_error("query_type should be either UCSC or BED");
  }

  const hictk::balancing::Method normalization_method{normalization.value_or("NONE")};

  if (normalization_method != hictk::balancing::Method::NONE()) {
    count_type = "float";
  }

  return fetch_impl(
      *f, range1, range2, normalization_method,
      std::visit([](const auto &ct) { return map_py_numeric_to_cpp_type(ct); }, count_type), join,
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED,
      diagonal_band_width);
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

static nb::dict get_attributes(const File &f) {
  if (f->is_cooler()) {
    return get_cooler_attrs(f->get<hictk::cooler::File>());
  }
  return get_hic_attrs(f->get<hictk::hic::File>());
}

static std::vector<std::string> avail_normalizations(const File &f) {
  const auto norms_ = f->avail_normalizations();
  std::vector<std::string> norms{norms_.size()};
  std::transform(norms_.begin(), norms_.end(), norms.begin(),
                 [](const auto &norm) { return norm.to_string(); });

  return norms;
}

// macro used to reduce boilerplate required to call methods on File objects
#define HICTKPY_CALL_METHOD_CHECKED(method) [](const File &x) { return x->method(); }

static auto get_chromosomes(const File &f, bool include_ALL) {
  return get_chromosomes_from_object(*f, include_ALL);
}

static auto get_bins(const File &f) { return get_bins_from_object(*f); }

static std::int64_t get_nchroms(const File &f, bool include_ALL) {
  return static_cast<std::int64_t>(f->nchroms(include_ALL));
}

static std::filesystem::path get_path(const File &f) { return std::filesystem::path{f->path()}; }

static auto get_weights(const File &f, std::string_view normalization, bool divisive) {
  using WeightVector = nb::ndarray<nb::numpy, nb::ndim<1>, nb::c_contig, double>;
  if (normalization == "NONE") {
    return std::optional<WeightVector>{};
  }

  const auto type = divisive ? hictk::balancing::Weights::Type::DIVISIVE
                             : hictk::balancing::Weights::Type::MULTIPLICATIVE;

  return std::make_optional(make_owning_numpy(f->normalization(normalization).to_vector(type)));
}

static nb::object get_weights_df(const File &f, const std::vector<std::string> &normalizations,
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
        f->nbins(), arrow::Buffer::FromVector(f->normalization(normalization).to_vector(type)),
        nullptr, 0, 0));
  }

  auto schema = std::make_shared<arrow::Schema>(std::move(fields));

  return export_pyarrow_table(arrow::Table::Make(std::move(schema), columns))
      .attr("to_pandas")(nb::arg("self_destruct") = true);
}

static bool has_normalization(const File &f, std::string_view name) {
  return f->has_normalization(name);
}

static std::string repr(const File &f) { return fmt::format(FMT_STRING("File({})"), f->uri()); }

static File &ctx_enter(File &f) { return f; }

static void ctx_exit(File &f, [[maybe_unused]] nb::handle exc_type,
                     [[maybe_unused]] nb::handle exc_value, [[maybe_unused]] nb::handle traceback) {
  f.try_close();
}

[[noreturn]] static void throw_closed_file_exc(std::string_view uri) {
  throw std::runtime_error(fmt::format(
      FMT_STRING("caught an attempt to access file \"{}\", which has already been closed"), uri));
}

[[nodiscard]] static std::optional<std::uint32_t> sanitize_resolution(
    std::optional<std::int32_t> resolution) {
  if (!resolution.has_value()) {
    return {};
  }

  if (resolution < 0) {
    throw std::invalid_argument("resolution must be non-negative");
  }

  return static_cast<std::uint32_t>(*resolution);
}

[[nodiscard]] static std::string get_uri_ts(const hictk::File &f) {
  if (f.is_cooler()) {
    HICTKPY_LOCK_COOLER_MTX_SCOPED
    return f.uri();
  }
  return f.uri();
}

File::File(hictk::File f) : _fp(std::move(f)), _uri(get_uri_ts(*_fp)) {}

File::File(hictk::cooler::File f) : File(hictk::File(std::move(f))) {}

File::File(hictk::hic::File f) : File(hictk::File(std::move(f))) {}

[[nodiscard]] static hictk::File open_file_ts(const std::filesystem::path &path,
                                              std::optional<std::int32_t> resolution,
                                              std::string_view matrix_type,
                                              std::string_view matrix_unit) {
  HICTKPY_LOCK_COOLER_MTX_SCOPED
  return hictk::File{path.string(), sanitize_resolution(resolution),
                     hictk::hic::ParseMatrixTypeStr(std::string{matrix_type}),
                     hictk::hic::ParseUnitStr(std::string{matrix_unit})};
}

File::File(const std::filesystem::path &path, std::optional<std::int32_t> resolution,
           std::string_view matrix_type, std::string_view matrix_unit)
    : File(open_file_ts(path, resolution, matrix_type, matrix_unit)) {}

File::~File() noexcept { std::ignore = try_close(); }

hictk::File *File::operator->() { return &**this; }

const hictk::File *File::operator->() const { return &**this; }

hictk::File &File::operator*() {
  if (_fp.has_value()) {
    return *_fp;
  }
  throw_closed_file_exc(_uri);
}

const hictk::File &File::operator*() const {
  if (_fp.has_value()) {
    return *_fp;
  }
  throw_closed_file_exc(_uri);
}

void File::close() {
  if (_fp.has_value()) {
    [[maybe_unused]] const auto lck = lock();
    _fp.reset();
  }
}

bool File::try_close() noexcept {
  if (!_fp) {
    return true;
  }

  try {
    try {
      close();
      return true;
    } catch (const std::exception &e) {
      raise_python_runtime_warning(FMT_STRING("an error occurred while closing file \"{}\": {}"),
                                   _uri, e.what());
    } catch (...) {
      raise_python_runtime_warning(
          FMT_STRING("an error occurred while closing file \"{}\": unknown error"), _uri);
    }
  } catch (...) {  // NOLINT
  }

  return false;
}

bool File::is_cooler(const std::filesystem::path &uri) {
  HICTKPY_LOCK_COOLER_MTX_SCOPED
  return static_cast<bool>(hictk::cooler::utils::is_cooler(uri.string()));
}

bool File::is_hic(const std::filesystem::path &uri) { return hictk::hic::utils::is_hic_file(uri); }

CoolerGlobalLock::UniqueLock File::lock() const {
  if (_fp.has_value() && _fp->is_cooler()) {
    return CoolerGlobalLock::lock();
  }
  return {};
}

void File::bind(nb::module_ &m) {
  auto file =
      nb::class_<File>(m, "File", "Class representing a file handle to a .cool or .hic file.");

  file.def(nb::init<const std::filesystem::path &, std::optional<std::int32_t>, std::string_view,
                    std::string_view>(),
           nb::arg("path"), nb::arg("resolution") = nb::none(), nb::arg("matrix_type") = "observed",
           nb::arg("matrix_unit") = "BP",
           "Construct a file object to a .hic, .cool or .mcool file given the file path and "
           "resolution.\n"
           "Resolution is ignored when opening single-resolution Cooler files.");

  file.def("__repr__", &repr, nb::rv_policy::move);

  file.def("__enter__", &ctx_enter, nb::rv_policy::reference_internal);

  file.def("__exit__", &ctx_exit,
           // clang-format off
           nb::arg("exc_type") = nb::none(),
           nb::arg("exc_value") = nb::none(),
           nb::arg("traceback") = nb::none()
           // clang-format on
  );

  file.def("uri", HICTKPY_CALL_METHOD_CHECKED(uri), "Return the file URI.", nb::rv_policy::move);
  file.def("path", &get_path, "Return the file path.", nb::rv_policy::move);

  file.def("is_hic", HICTKPY_CALL_METHOD_CHECKED(is_hic), "Test whether file is in .hic format.");
  file.def("is_cooler", HICTKPY_CALL_METHOD_CHECKED(is_cooler),
           "Test whether file is in .cool format.");

  file.def("close", &File::close, "Manually close the file handle.");

  file.def("chromosomes", &get_chromosomes, nb::arg("include_ALL") = false,
           "Get chromosome sizes as a dictionary mapping names to sizes.",
           nb::rv_policy::take_ownership);
  file.def("bins", &get_bins, "Get table of bins.", nb::sig("def bins(self) -> hictkpy.BinTable"),
           nb::rv_policy::move);

  file.def("resolution", HICTKPY_CALL_METHOD_CHECKED(resolution), "Get the bin size in bp.");
  file.def("nbins", HICTKPY_CALL_METHOD_CHECKED(nbins), "Get the total number of bins.");
  file.def("nchroms", &get_nchroms, nb::arg("include_ALL") = false,
           "Get the total number of chromosomes.");

  file.def("attributes", &get_attributes, "Get file attributes as a dictionary.",
           nb::rv_policy::take_ownership);

  file.def("fetch", &fetch, nb::keep_alive<0, 1>(), nb::arg("range1") = nb::none(),
           nb::arg("range2") = nb::none(), nb::arg("normalization") = nb::none(),
           nb::arg("count_type") = "int32", nb::arg("join") = false, nb::arg("query_type") = "UCSC",
           nb::arg("diagonal_band_width") = nb::none(),
           "Fetch interactions overlapping a region of interest.", nb::rv_policy::move);

  file.def("avail_normalizations", &avail_normalizations,
           "Get the list of available normalizations.", nb::rv_policy::move);
  file.def("has_normalization", &has_normalization, nb::arg("normalization"),
           "Check whether a given normalization is available.");
  file.def("weights", &get_weights, nb::arg("name"), nb::arg("divisive") = true,
           "Fetch the balancing weights for the given normalization method.",
           nb::rv_policy::take_ownership);
  file.def(
      "weights", &get_weights_df, nb::arg("names"), nb::arg("divisive") = true,
      "Fetch the balancing weights for the given normalization methods."
      "Weights are returned as a pandas.DataFrame.",
      nb::sig("def weights(self, names: collections.abc.Sequence[str], divisive: bool = True) -> "
              "pandas.DataFrame"),
      nb::rv_policy::take_ownership);
}

void cooler::init_global_state() {
  // this should be called early on during module declaration, and it is mostly useful to avoid
  // false positives from the sanitizers
  HICTKPY_LOCK_COOLER_MTX_SCOPED
  if (H5open() < 0) {
    throw std::runtime_error("failed to initialize HDF5 library!");
  }
}

}  // namespace hictkpy
