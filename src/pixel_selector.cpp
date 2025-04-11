// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef _WIN32
// Workaround bug several symbol redefinition errors due to something including <winsock.h>
#include <winsock2.h>
#endif

#include <arrow/table.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/cooler/pixel_selector.hpp>
#include <hictk/fmt.hpp>
#include <hictk/hic/pixel_selector.hpp>
#include <hictk/pixel.hpp>
#include <hictk/transformers/common.hpp>
#include <hictk/transformers/diagonal_band.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <hictk/transformers/to_dataframe.hpp>
#include <hictk/transformers/to_dense_matrix.hpp>
#include <hictk/transformers/to_sparse_matrix.hpp>
#include <hictkpy/common.hpp>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel_aggregator.hpp"
#include "hictkpy/pixel_selector.hpp"
#include "hictkpy/to_pyarrow.hpp"

namespace nb = nanobind;

namespace hictkpy {

static std::optional<std::uint64_t> transform_diagonal_band_width(std::optional<std::int64_t> w) {
  if (!w.has_value()) {
    return {};
  }

  if (*w < 0) {
    throw std::invalid_argument("diagonal_band_width cannot be negative");
  }

  return static_cast<std::uint64_t>(*w);
}

PixelSelector::PixelSelector(std::shared_ptr<const hictk::cooler::PixelSelector> sel_,
                             std::string_view type, bool join,
                             std::optional<std::int64_t> diagonal_band_width)
    : selector(std::move(sel_)),
      pixel_count(parse_count_type(type)),
      pixel_format(join ? PixelFormat::BG2 : PixelFormat::COO),
      _diagonal_band_width(transform_diagonal_band_width(diagonal_band_width)) {}

PixelSelector::PixelSelector(std::shared_ptr<const hictk::hic::PixelSelector> sel_,
                             std::string_view type, bool join,
                             std::optional<std::int64_t> diagonal_band_width)
    : selector(std::move(sel_)),
      pixel_count(parse_count_type(type)),
      pixel_format(join ? PixelFormat::BG2 : PixelFormat::COO),
      _diagonal_band_width(transform_diagonal_band_width(diagonal_band_width)) {}

PixelSelector::PixelSelector(std::shared_ptr<const hictk::hic::PixelSelectorAll> sel_,
                             std::string_view type, bool join,
                             std::optional<std::int64_t> diagonal_band_width)
    : selector(std::move(sel_)),
      pixel_count(parse_count_type(type)),
      pixel_format(join ? PixelFormat::BG2 : PixelFormat::COO),
      _diagonal_band_width(transform_diagonal_band_width(diagonal_band_width)) {}

std::string PixelSelector::repr() const {
  if (!coord1()) {
    return fmt::format(FMT_STRING("PixelSelector(ALL; {}; {})"),
                       pixel_format == PixelFormat::COO ? "COO" : "BG2",
                       count_type_to_str(pixel_count));
  }

  return fmt::format(FMT_STRING("PixelSelector({}:{}-{}; {}:{}-{}; {}; {})"),
                     coord1().bin1.chrom().name(), coord1().bin1.start(), coord1().bin2.end(),
                     coord2().bin1.chrom().name(), coord2().bin1.start(), coord2().bin2.end(),
                     pixel_format == PixelFormat::COO ? "COO" : "BG2",
                     count_type_to_str(pixel_count));
}

// NOLINTNEXTLINE(bugprone-exception-escape)
hictk::PixelCoordinates PixelSelector::coord1() const noexcept {
  assert(!selector.valueless_by_exception());
  return std::visit(
      [](const auto& s) -> hictk::PixelCoordinates {
        if constexpr (std::is_same_v<std::decay_t<decltype(*s)>, hictk::hic::PixelSelectorAll>) {
          return {};
        } else {
          return s->coord1();
        }
      },
      selector);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
hictk::PixelCoordinates PixelSelector::coord2() const noexcept {
  assert(!selector.valueless_by_exception());
  return std::visit(
      [](const auto& s) -> hictk::PixelCoordinates {
        if constexpr (std::is_same_v<std::decay_t<decltype(*s)>, hictk::hic::PixelSelectorAll>) {
          return {};
        } else {
          return s->coord2();
        }
      },
      selector);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
const hictk::BinTable& PixelSelector::bins() const noexcept {
  assert(!selector.valueless_by_exception());
  return std::visit([](const auto& s) -> const hictk::BinTable& { return s->bins(); }, selector);
}

[[nodiscard]] static PixelSelector::GenomicCoordTuple coords_to_tuple(
    const hictk::PixelCoordinates& coords, const hictk::BinTable& bins) {
  if (!coords) {
    return {"ALL", 0, static_cast<std::int64_t>(bins.size())};
  }

  assert(coords.bin1.chrom() == coords.bin2.chrom());

  return {std::string{coords.bin1.chrom().name()}, static_cast<std::int64_t>(coords.bin1.start()),
          static_cast<std::int64_t>(coords.bin2.end())};
}

auto PixelSelector::get_coord1() const -> GenomicCoordTuple {
  return coords_to_tuple(coord1(), bins());
}

auto PixelSelector::get_coord2() const -> GenomicCoordTuple {
  return coords_to_tuple(coord2(), bins());
}

std::int64_t PixelSelector::size(bool upper_triangular) const {
  return std::visit(
      [&](const auto& sel_ptr) {
        return static_cast<std::int64_t>(sel_ptr->size(upper_triangular));
      },
      selector);
}

template <typename N, typename PixelSelector>
[[nodiscard]] static nb::iterator make_bg2_iterable(
    const PixelSelector& sel, std::optional<std::uint64_t> diagonal_band_width) {
  if constexpr (std::is_floating_point_v<N> && !std::is_same_v<N, double>) {
    return make_bg2_iterable<double>(sel, diagonal_band_width);
  } else if constexpr (std::is_integral_v<N> && !std::is_same_v<N, std::int64_t>) {
    return make_bg2_iterable<std::int64_t>(sel, diagonal_band_width);
  }

  if (diagonal_band_width.has_value()) {
    const hictk::transformers::DiagonalBand band_sel(sel.template begin<N>(), sel.template end<N>(),
                                                     *diagonal_band_width);
    const hictk::transformers::JoinGenomicCoords jsel{band_sel.begin(), band_sel.end(),
                                                      sel.bins_ptr()};
    return nb::make_iterator(nb::type<hictkpy::PixelSelector>(), "PixelIterator", jsel.begin(),
                             jsel.end());
  }

  const hictk::transformers::JoinGenomicCoords jsel{sel.template begin<N>(), sel.template end<N>(),
                                                    sel.bins_ptr()};
  return nb::make_iterator(nb::type<hictkpy::PixelSelector>(), "PixelIterator", jsel.begin(),
                           jsel.end());
}

template <typename N, typename PixelSelector>
[[nodiscard]] static nb::iterator make_coo_iterable(
    const PixelSelector& sel, std::optional<std::uint64_t> diagonal_band_width) {
  if constexpr (std::is_floating_point_v<N> && !std::is_same_v<N, double>) {
    return make_coo_iterable<double>(sel, diagonal_band_width);
  } else if constexpr (std::is_integral_v<N> && !std::is_same_v<N, std::int64_t>) {
    return make_coo_iterable<std::int64_t>(sel, diagonal_band_width);
  }

  if (diagonal_band_width.has_value()) {
    const hictk::transformers::DiagonalBand band_sel(sel.template begin<N>(), sel.template end<N>(),
                                                     *diagonal_band_width);
    return nb::make_iterator(nb::type<hictkpy::PixelSelector>(), "PixelIterator", band_sel.begin(),
                             band_sel.end());
  }

  return nb::make_iterator(nb::type<hictkpy::PixelSelector>(), "PixelIterator",
                           sel.template begin<N>(), sel.template end<N>());
}

nb::iterator PixelSelector::make_iterable() const {
  return std::visit(
      [&](const auto& sel_ptr) -> nb::iterator {
        assert(!!sel_ptr);
        return std::visit(
            [&]([[maybe_unused]] auto count) -> nb::iterator {
              using N = decltype(count);
              if (pixel_format == PixelFormat::BG2) {
                return make_bg2_iterable<N>(*sel_ptr, _diagonal_band_width);
              }
              assert(pixel_format == PixelFormat::COO);
              return make_coo_iterable<N>(*sel_ptr, _diagonal_band_width);
            },
            pixel_count);
      },
      selector);
}

template <typename N, typename PixelSelector>
[[nodiscard]] static std::shared_ptr<arrow::Table> make_bg2_arrow_df(
    const PixelSelector& sel, hictk::transformers::QuerySpan span,
    std::optional<std::uint64_t> diagonal_band_width) {
  if constexpr (std::is_same_v<N, long double>) {
    return make_bg2_arrow_df<double>(sel, span, diagonal_band_width);
  } else {
    return hictk::transformers::ToDataFrame(
        sel.template begin<N>(), sel.template end<N>(), hictk::transformers::DataFrameFormat::BG2,
        sel.bins_ptr(), span, false, true, 256'000, diagonal_band_width)();
  }
}

template <typename N, typename PixelSelector>
[[nodiscard]] static std::shared_ptr<arrow::Table> make_coo_arrow_df(
    const PixelSelector& sel, hictk::transformers::QuerySpan span,
    std::optional<std::uint64_t> diagonal_band_width) {
  if constexpr (std::is_same_v<N, long double>) {
    return make_coo_arrow_df<double>(sel, span, diagonal_band_width);
  } else {
    return hictk::transformers::ToDataFrame(
        sel.template begin<N>(), sel.template end<N>(), hictk::transformers::DataFrameFormat::COO,
        sel.bins_ptr(), span, false, true, 256'000, diagonal_band_width)();
  }
}

nb::object PixelSelector::to_arrow(std::string_view span) const {
  std::ignore = import_pyarrow_checked();

  const auto query_span = parse_span(span);
  auto table = std::visit(
      [&](const auto& sel_ptr) -> std::shared_ptr<arrow::Table> {
        assert(!!sel_ptr);
        return std::visit(
            [&]([[maybe_unused]] auto count) -> std::shared_ptr<arrow::Table> {
              using N = decltype(count);
              if (pixel_format == PixelFormat::BG2) {
                return make_bg2_arrow_df<N>(*sel_ptr, query_span, _diagonal_band_width);
              }
              assert(pixel_format == PixelFormat::COO);
              return make_coo_arrow_df<N>(*sel_ptr, query_span, _diagonal_band_width);
            },
            pixel_count);
      },
      selector);

  return export_pyarrow_table(std::move(table));
}

nb::object PixelSelector::to_pandas(std::string_view span) const {
  import_module_checked("pandas");
  return to_arrow(span).attr("to_pandas")(nb::arg("self_destruct") = true);
}

nb::object PixelSelector::to_df(std::string_view span) const { return to_pandas(span); }

template <typename N, typename PixelSelector>
[[nodiscard]] static nb::object make_csr_matrix(std::shared_ptr<const PixelSelector> sel,
                                                hictk::transformers::QuerySpan span,
                                                std::optional<std::uint64_t> diagonal_band_width) {
  if constexpr (std::is_same_v<N, long double>) {
    return make_csr_matrix<double>(std::move(sel), span, diagonal_band_width);
  } else {
    return nb::cast(hictk::transformers::ToSparseMatrix(std::move(sel), N{}, span, false,
                                                        diagonal_band_width)());
  }
}

nb::object PixelSelector::to_csr(std::string_view span) const {
  import_module_checked("scipy");
  const auto query_span = parse_span(span);

  return std::visit(
      [&](auto sel_ptr) -> nb::object {
        return std::visit(
            [&]([[maybe_unused]] auto count) -> nb::object {
              using N = decltype(count);
              return make_csr_matrix<N>(std::move(sel_ptr), query_span, _diagonal_band_width);
            },
            pixel_count);
      },
      selector);
}

nb::object PixelSelector::to_coo(std::string_view span) const {
  import_module_checked("scipy");
  return to_csr(span).attr("tocoo")(false);
}

template <typename N, typename PixelSelector>
[[nodiscard]] static nb::object make_numpy_matrix(
    std::shared_ptr<const PixelSelector> sel, hictk::transformers::QuerySpan span,
    std::optional<std::uint64_t> diagonal_band_width) {
  if constexpr (std::is_same_v<N, long double>) {
    return make_numpy_matrix<double>(std::move(sel), span, diagonal_band_width);
  } else {
    return nb::cast(
        hictk::transformers::ToDenseMatrix(std::move(sel), N{}, span, diagonal_band_width)());
  }
}

nb::object PixelSelector::to_numpy(std::string_view span) const {
  std::ignore = import_module_checked("numpy");

  const auto query_span = parse_span(span);

  return std::visit(
      [&](auto sel_ptr) -> nb::object {
        return std::visit(
            [&]([[maybe_unused]] auto count) -> nb::object {
              using N = decltype(count);
              return make_numpy_matrix<N>(std::move(sel_ptr), query_span, _diagonal_band_width);
            },
            pixel_count);
      },
      selector);
}

template <typename PixelIt>
[[nodiscard]] static Stats aggregate_pixels(PixelIt first, PixelIt last, std::uint64_t size,
                                            bool keep_nans, bool keep_infs, bool keep_zeros,
                                            bool exact,
                                            const phmap::flat_hash_set<std::string>& metrics) {
  static_assert(!std::is_same_v<PixelSelector, hictk::PixelSelector>);
  if (keep_nans && keep_infs) {
    return PixelAggregator<PixelIt>{}.template compute<true, true>(
        std::move(first), std::move(last), size, metrics, keep_zeros, exact);
  }
  if (keep_nans) {
    return PixelAggregator<PixelIt>{}.template compute<true, false>(
        std::move(first), std::move(last), size, metrics, keep_zeros, exact);
  }
  if (keep_infs) {
    return PixelAggregator<PixelIt>{}.template compute<false, true>(
        std::move(first), std::move(last), size, metrics, keep_zeros, exact);
  }
  return PixelAggregator<PixelIt>{}.template compute<false, false>(
      std::move(first), std::move(last), size, metrics, keep_zeros, exact);
}

[[nodiscard]] static Stats aggregate_pixels(const PixelSelector::SelectorVar& sel,
                                            const PixelSelector::PixelVar& count, bool keep_nans,
                                            bool keep_infs, bool keep_zeros, bool exact,
                                            std::optional<std::uint64_t> diagonal_band_width,
                                            const phmap::flat_hash_set<std::string>& metrics) {
  // MSVC gets confused if fixed_bin_size is declared inside the nested call to std::visit down
  // below
  const auto fixed_bin_size = std::visit(
      [&](const auto& sel_ptr) {
        assert(sel_ptr);
        return sel_ptr->bins().type() == hictk::BinTable::Type::fixed;
      },
      sel);

  if (!fixed_bin_size && keep_zeros) {
    throw std::runtime_error(
        "calculating statistics including zeros on files with bin tables other than "
        "\"fixed\" bin size is not supported.");
  }

  if (diagonal_band_width.has_value() && keep_zeros) {
    // TODO look into how difficult it is to support this
    throw std::runtime_error(
        "using diagonal_band_width is not currently supported when keep_zeros=True");
  }

  // All the explicit captures are required to make MSVC happy
  return std::visit(
      [&metrics, &count, fixed_bin_size, keep_nans, keep_infs, keep_zeros, diagonal_band_width,
       exact](const auto& sel_ptr) {
        assert(sel_ptr);
        return std::visit(
            [&metrics, &sel_ptr, fixed_bin_size, keep_nans, keep_infs, keep_zeros,
             diagonal_band_width, exact]([[maybe_unused]] const auto& count_) {
              using N = remove_cvref_t<decltype(count_)>;

              auto first = sel_ptr->template begin<N>();
              auto last = sel_ptr->template end<N>();

              if (diagonal_band_width.has_value()) {
                hictk::transformers::DiagonalBand sel_diag{
                    std::move(first), std::move(last),
                    static_cast<std::uint64_t>(*diagonal_band_width)};
                return aggregate_pixels(sel_diag.begin(), sel_diag.end(),
                                        fixed_bin_size ? sel_ptr->size() : 0, keep_nans, keep_infs,
                                        keep_zeros, exact, metrics);
              }
              return aggregate_pixels(std::move(first), std::move(last),
                                      fixed_bin_size ? sel_ptr->size() : 0, keep_nans, keep_infs,
                                      keep_zeros, exact, metrics);
            },
            count);
      },
      sel);
}

nb::dict PixelSelector::describe(const std::vector<std::string>& metrics, bool keep_nans,
                                 bool keep_infs, bool keep_zeros, bool exact) const {
  const auto stats =
      aggregate_pixels(selector, pixel_count, keep_nans, keep_infs, keep_zeros, exact,
                       _diagonal_band_width, {metrics.begin(), metrics.end()});

  using StatsDict =
      nanobind::typed<nanobind::dict, std::string, std::variant<std::int64_t, double>>;

  StatsDict stats_py{};

  for (const auto& metric : metrics) {
    stats_py[nb::cast(metric)] = nb::none();
  }

  if (stats.nnz) {
    stats_py["nnz"] = *stats.nnz;
  }
  if (stats.sum) {
    stats_py["sum"] =
        std::visit([](const auto n) -> nb::object { return nb::cast(n); }, *stats.sum);
  }
  if (stats.min) {
    stats_py["min"] =
        std::visit([](const auto n) -> nb::object { return nb::cast(n); }, *stats.min);
  }
  if (stats.max) {
    stats_py["max"] =
        std::visit([](const auto n) -> nb::object { return nb::cast(n); }, *stats.max);
  }
  if (stats.mean) {
    stats_py["mean"] = *stats.mean;
  }
  if (stats.variance) {
    stats_py["variance"] = *stats.variance;
  }
  if (stats.skewness) {
    stats_py["skewness"] = *stats.skewness;
  }
  if (stats.kurtosis) {
    stats_py["kurtosis"] = *stats.kurtosis;
  }

  return stats_py;
}

std::int64_t PixelSelector::nnz(bool keep_nans, bool keep_infs) const {
  return *aggregate_pixels(selector, pixel_count, keep_nans, keep_infs, false, false,
                           _diagonal_band_width, {"nnz"})
              .nnz;
}

nb::object PixelSelector::sum(bool keep_nans, bool keep_infs) const {
  const auto stats = aggregate_pixels(selector, pixel_count, keep_nans, keep_infs, false, false,
                                      _diagonal_band_width, {"sum"});
  return std::visit([](const auto n) -> nb::object { return nb::cast(n); }, *stats.sum);
}

nb::object PixelSelector::min(bool keep_nans, bool keep_infs, bool keep_zeros) const {
  const auto stats = aggregate_pixels(selector, pixel_count, keep_nans, keep_infs, keep_zeros,
                                      false, _diagonal_band_width, {"min"});
  return std::visit([](const auto n) -> nb::object { return nb::cast(n); }, *stats.min);
}

nb::object PixelSelector::max(bool keep_nans, bool keep_infs, bool keep_zeros) const {
  const auto stats = aggregate_pixels(selector, pixel_count, keep_nans, keep_infs, keep_zeros,
                                      false, _diagonal_band_width, {"max"});
  return std::visit([](const auto n) -> nb::object { return nb::cast(n); }, *stats.max);
}

double PixelSelector::mean(bool keep_nans, bool keep_infs, bool keep_zeros) const {
  return *aggregate_pixels(selector, pixel_count, keep_nans, keep_infs, keep_zeros, false,
                           _diagonal_band_width, {"mean"})
              .mean;
}

double PixelSelector::variance(bool keep_nans, bool keep_infs, bool keep_zeros, bool exact) const {
  return *aggregate_pixels(selector, pixel_count, keep_nans, keep_infs, keep_zeros, exact,
                           _diagonal_band_width, {"variance"})
              .variance;
}

double PixelSelector::skewness(bool keep_nans, bool keep_infs, bool keep_zeros, bool exact) const {
  return *aggregate_pixels(selector, pixel_count, keep_nans, keep_infs, keep_zeros, exact,
                           _diagonal_band_width, {"skewness"})
              .skewness;
}

double PixelSelector::kurtosis(bool keep_nans, bool keep_infs, bool keep_zeros, bool exact) const {
  return *aggregate_pixels(selector, pixel_count, keep_nans, keep_infs, keep_zeros, exact,
                           _diagonal_band_width, {"kurtosis"})
              .kurtosis;
}

auto PixelSelector::parse_span(std::string_view span) -> QuerySpan {
  if (span == "upper_triangle") {
    return QuerySpan::upper_triangle;
  }
  if (span == "lower_triangle") {
    return QuerySpan::lower_triangle;
  }
  if (span == "full") {
    return QuerySpan::full;
  }

  throw std::runtime_error(
      fmt::format(FMT_STRING("unrecognized query span \"{}\". Supported query spans are: "
                             "upper_triangle, lower_triangle, and full"),
                  span));
}

auto PixelSelector::parse_count_type(std::string_view type) -> PixelVar {
  return map_dtype_to_type(type);
}

std::string_view PixelSelector::count_type_to_str(const PixelVar& var) {
  // NOLINTBEGIN(*-avoid-magic-numbers)
  static_assert(sizeof(float) == 4);
  static_assert(sizeof(double) == 8);
  // NOLINTEND(*-avoid-magic-numbers)

  if (std::holds_alternative<std::uint8_t>(var)) {
    return "uint8";
  }
  if (std::holds_alternative<std::uint16_t>(var)) {
    return "uint16";
  }
  if (std::holds_alternative<std::uint32_t>(var)) {
    return "uint32";
  }
  if (std::holds_alternative<std::uint64_t>(var)) {
    return "uint64";
  }
  if (std::holds_alternative<std::int8_t>(var)) {
    return "int8";
  }
  if (std::holds_alternative<std::int16_t>(var)) {
    return "int16";
  }
  if (std::holds_alternative<std::int32_t>(var)) {
    return "int32";
  }
  if (std::holds_alternative<std::int64_t>(var)) {
    return "int64";
  }
  if (std::holds_alternative<float>(var)) {
    return "float32";
  }
  if (std::holds_alternative<double>(var)) {
    return "float64";
  }

  unreachable_code();
}

void PixelSelector::bind(nb::module_& m) {
  auto sel = nb::class_<PixelSelector>(
      m, "PixelSelector",
      "Class representing pixels overlapping with the given genomic intervals.");

  sel.def(nb::init<std::shared_ptr<const hictk::cooler::PixelSelector>, std::string_view, bool,
                   std::optional<std::int64_t>>(),
          nb::arg("selector"), nb::arg("type"), nb::arg("join"),
          nb::arg("diagonal_band_width") = nb::none());
  sel.def(nb::init<std::shared_ptr<const hictk::hic::PixelSelector>, std::string_view, bool,
                   std::optional<std::int64_t>>(),
          nb::arg("selector"), nb::arg("type"), nb::arg("join"),
          nb::arg("diagonal_band_width") = nb::none());
  sel.def(nb::init<std::shared_ptr<const hictk::hic::PixelSelectorAll>, std::string_view, bool,
                   std::optional<std::int64_t>>(),
          nb::arg("selector"), nb::arg("type"), nb::arg("join"),
          nb::arg("diagonal_band_width") = nb::none());

  sel.def("__repr__", &PixelSelector::repr, nb::rv_policy::move);

  sel.def("coord1", &PixelSelector::get_coord1, "Get query coordinates for the first dimension.",
          nb::rv_policy::move);
  sel.def("coord2", &PixelSelector::get_coord2, "Get query coordinates for the second dimension.",
          nb::rv_policy::move);
  sel.def("size", &PixelSelector::size, nb::arg("upper_triangular") = true,
          "Get the number of pixels overlapping with the given query.");

  sel.def("__iter__", &PixelSelector::make_iterable, nb::keep_alive<0, 1>(),
          nb::sig("def __iter__(self) -> hictkpy.PixelIterator"),
          "Return an iterator over the selected pixels.", nb::rv_policy::take_ownership);

  sel.def("to_arrow", &PixelSelector::to_arrow, nb::arg("query_span") = "upper_triangle",
          nb::sig("def to_arrow(self, query_span: str = \"upper_triangle\") -> pyarrow.Table"),
          "Retrieve interactions as a pyarrow.Table.", nb::rv_policy::take_ownership);
  sel.def("to_pandas", &PixelSelector::to_pandas, nb::arg("query_span") = "upper_triangle",
          nb::sig("def to_pandas(self, query_span: str = \"upper_triangle\") -> pandas.DataFrame"),
          "Retrieve interactions as a pandas DataFrame.", nb::rv_policy::take_ownership);
  sel.def("to_df", &PixelSelector::to_df, nb::arg("query_span") = "upper_triangle",
          nb::sig("def to_df(self, query_span: str = \"upper_triangle\") -> pandas.DataFrame"),
          "Alias to to_pandas().", nb::rv_policy::take_ownership);
  sel.def("to_numpy", &PixelSelector::to_numpy, nb::arg("query_span") = "full",
          nb::sig("def to_numpy(self, query_span: str = \"full\") -> numpy.ndarray"),
          "Retrieve interactions as a numpy 2D matrix.", nb::rv_policy::move);
  sel.def(
      "to_coo", &PixelSelector::to_coo, nb::arg("query_span") = "upper_triangle",
      nb::sig("def to_coo(self, query_span: str = \"upper_triangle\") -> scipy.sparse.coo_matrix"),
      "Retrieve interactions as a SciPy COO matrix.", nb::rv_policy::take_ownership);
  sel.def(
      "to_csr", &PixelSelector::to_csr, nb::arg("query_span") = "upper_triangle",
      nb::sig("def to_csr(self, query_span: str = \"upper_triangle\") -> scipy.sparse.csr_matrix"),
      "Retrieve interactions as a SciPy CSR matrix.", nb::rv_policy::move);

  using PixelIt = hictk::ThinPixel<int>*;
  static const std::vector<std::string> known_metrics(
      PixelAggregator<PixelIt>::valid_metrics.begin(),
      PixelAggregator<PixelIt>::valid_metrics.end());
  static const auto describe_cmd_help = fmt::format(
      FMT_STRING(
          "Compute one or more descriptive metrics in the most efficient way possible. "
          "Known metrics: {}. "
          "When a metric cannot be computed (e.g. because metrics=[\"variance\"], "
          "but selector overlaps with a single pixel), the value for that metric is set to None. "
          "When keep_infs or keep_nans are set to True, and keep_zeros=True, nan and/or inf "
          "values are treated as zeros. "
          "By default, metrics are estimated by doing a single pass through the data. "
          "The estimates are stable and usually very accurate. "
          "However, if you require exact values, you can specify exact=True."),
      fmt::join(known_metrics, ", "));
  sel.def("describe", &PixelSelector::describe, nb::arg("metrics") = known_metrics,
          nb::arg("keep_nans") = false, nb::arg("keep_infs") = false, nb::arg("keep_zeros") = false,
          nb::arg("exact") = false, describe_cmd_help.c_str());
  sel.def("nnz", &PixelSelector::nnz, nb::arg("keep_nans") = false, nb::arg("keep_infs") = false,
          "Get the number of non-zero entries for the current pixel selection. See documentation "
          "for describe() for more details.");
  sel.def("sum", &PixelSelector::sum, nb::arg("keep_nans") = false, nb::arg("keep_infs") = false,
          nb::sig("def sum(self, keep_nans: bool = False, keep_infs: bool = False) -> int | float"),
          "Get the total number of interactions for the current pixel selection. See documentation "
          "for describe() for more details.");
  sel.def(
      "min", &PixelSelector::min, nb::arg("keep_nans") = false, nb::arg("keep_infs") = false,
      nb::arg("keep_zeros") = false,
      nb::sig("def min(self, keep_nans: bool = False, keep_infs: bool = False, keep_zeros: bool = "
              "False) -> int | float | None"),
      "Get the minimum number of interactions for the current pixel selection. See documentation "
      "for describe() for more details.");
  sel.def(
      "max", &PixelSelector::max, nb::arg("keep_nans") = false, nb::arg("keep_infs") = false,
      nb::arg("keep_zeros") = false,
      nb::sig("def max(self, keep_nans: bool = False, keep_infs: bool = False, keep_zeros: bool = "
              "False) -> int | float | None"),
      "Get the maximum number of interactions for the current pixel selection. See documentation "
      "for describe() for more details.");
  sel.def(
      "mean", &PixelSelector::mean, nb::arg("keep_nans") = false, nb::arg("keep_infs") = false,
      nb::arg("keep_zeros") = false,
      nb::sig("def mean(self, keep_nans: bool = False, keep_infs: bool = False, keep_zeros: bool = "
              "False) -> float | None"),
      "Get the average number of interactions for the current pixel selection. See documentation "
      "for describe() for more details.");
  sel.def(
      "variance", &PixelSelector::variance, nb::arg("keep_nans") = false,
      nb::arg("keep_infs") = false, nb::arg("keep_zeros") = false, nb::arg("exact") = false,
      nb::sig("def variance(self, keep_nans: bool = False, keep_infs: bool = False, keep_zeros: "
              "bool = False, exact: bool = False) -> float | None"),
      "Get the variance of the number of interactions for the current pixel selection. See "
      "documentation for describe() for more details.");
  sel.def(
      "skewness", &PixelSelector::skewness, nb::arg("keep_nans") = false,
      nb::arg("keep_infs") = false, nb::arg("keep_zeros") = false, nb::arg("exact") = false,
      nb::sig("def skewness(self, keep_nans: bool = False, keep_infs: bool = False, keep_zeros: "
              "bool = False, exact: bool = False) -> float | None"),
      "Get the skewness of the number of interactions for the current pixel selection. See "
      "documentation for describe() for more details.");
  sel.def(
      "kurtosis", &PixelSelector::kurtosis, nb::arg("keep_nans") = false,
      nb::arg("keep_infs") = false, nb::arg("keep_zeros") = false, nb::arg("exact") = false,
      nb::sig("def kurtosis(self, keep_nans: bool = False, keep_infs: bool = False, keep_zeros: "
              "bool = False, exact: bool = False) -> float | None"),
      "Get the kurtosis of the number of interactions for the current pixel selection. See "
      "documentation for describe() for more details.");
}

}  // namespace hictkpy
