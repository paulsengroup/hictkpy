// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef _WIN32
// Workaround bug several symbol redefinition errors due to something including <winsock.h>
#include <winsock2.h>
#endif

#include <arrow/python/api.h>
#include <arrow/table.h>
#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/cooler/pixel_selector.hpp>
#include <hictk/fmt.hpp>
#include <hictk/hic/pixel_selector.hpp>
#include <hictk/transformers/common.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <hictk/transformers/to_dataframe.hpp>
#include <hictk/transformers/to_dense_matrix.hpp>
#include <hictk/transformers/to_sparse_matrix.hpp>
#include <hictkpy/common.hpp>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel_selector.hpp"

namespace nb = nanobind;

namespace hictkpy {

PixelSelector::PixelSelector(std::shared_ptr<const hictk::cooler::PixelSelector> sel_,
                             std::string_view type, bool join)
    : selector(std::move(sel_)),
      pixel_count(parse_count_type(type)),
      pixel_format(join ? PixelFormat::BG2 : PixelFormat::COO) {}

PixelSelector::PixelSelector(std::shared_ptr<const hictk::hic::PixelSelector> sel_,
                             std::string_view type, bool join)
    : selector(std::move(sel_)),
      pixel_count(parse_count_type(type)),
      pixel_format(join ? PixelFormat::BG2 : PixelFormat::COO) {}

PixelSelector::PixelSelector(std::shared_ptr<const hictk::hic::PixelSelectorAll> sel_,
                             std::string_view type, bool join)
    : selector(std::move(sel_)),
      pixel_count(parse_count_type(type)),
      pixel_format(join ? PixelFormat::BG2 : PixelFormat::COO) {}

std::string PixelSelector::repr() const {
  if (!coord1()) {
    return fmt::format(FMT_STRING("PixelSelector(ALL; {}; {})"),
                       pixel_format == PixelFormat::COO ? "COO" : "BG2",
                       count_type_to_str(pixel_count));
  }

  return fmt::format(FMT_STRING("PixelSelector({}, {}; {}; {})"), coord1(), coord2(),
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

auto PixelSelector::get_coord1() const -> PixelCoordTuple {
  const auto c = coord1();
  return PixelCoordTuple{std::make_tuple(c.bin1.chrom().name(), c.bin1.start(), c.bin1.end(),
                                         c.bin2.chrom().name(), c.bin2.start(), c.bin2.end())};
}

auto PixelSelector::get_coord2() const -> PixelCoordTuple {
  const auto c = coord2();
  return PixelCoordTuple{std::make_tuple(c.bin1.chrom().name(), c.bin1.start(), c.bin1.end(),
                                         c.bin2.chrom().name(), c.bin2.start(), c.bin2.end())};
}

template <typename N, typename PixelSelector>
[[nodiscard]] static nb::iterator make_bg2_iterable(const PixelSelector& sel) {
  if constexpr (std::is_floating_point_v<N> && !std::is_same_v<N, double>) {
    return make_bg2_iterable<double>(sel);
  } else if constexpr (std::is_integral_v<N> && !std::is_same_v<N, std::int64_t>) {
    return make_bg2_iterable<std::int64_t>(sel);
  }
  const hictk::transformers::JoinGenomicCoords jsel{sel.template begin<N>(), sel.template end<N>(),
                                                    sel.bins_ptr()};
  return nb::make_iterator(nb::type<hictkpy::PixelSelector>(), "PixelIterator", jsel.begin(),
                           jsel.end());
}

template <typename N, typename PixelSelector>
[[nodiscard]] static nb::iterator make_coo_iterable(const PixelSelector& sel) {
  if constexpr (std::is_floating_point_v<N> && !std::is_same_v<N, double>) {
    return make_coo_iterable<double>(sel);
  } else if constexpr (std::is_integral_v<N> && !std::is_same_v<N, std::int64_t>) {
    return make_coo_iterable<std::int64_t>(sel);
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
                return make_bg2_iterable<N>(*sel_ptr);
              }
              assert(pixel_format == PixelFormat::COO);
              return make_coo_iterable<N>(*sel_ptr);
            },
            pixel_count);
      },
      selector);
}

template <typename N, typename PixelSelector>
[[nodiscard]] static std::shared_ptr<arrow::Table> make_bg2_arrow_df(
    const PixelSelector& sel, hictk::transformers::QuerySpan span) {
  if constexpr (std::is_same_v<N, long double>) {
    return make_bg2_arrow_df<double>(sel, span);
  } else {
    return hictk::transformers::ToDataFrame(sel.template begin<N>(), sel.template end<N>(),
                                            hictk::transformers::DataFrameFormat::BG2,
                                            sel.bins_ptr(), span)();
  }
}

template <typename N, typename PixelSelector>
[[nodiscard]] static std::shared_ptr<arrow::Table> make_coo_arrow_df(
    const PixelSelector& sel, hictk::transformers::QuerySpan span) {
  if constexpr (std::is_same_v<N, long double>) {
    return make_coo_arrow_df<double>(sel, span);
  } else {
    return hictk::transformers::ToDataFrame(sel.template begin<N>(), sel.template end<N>(),
                                            hictk::transformers::DataFrameFormat::COO,
                                            sel.bins_ptr(), span)();
  }
}

nb::object PixelSelector::to_arrow(std::string_view span) const {
  const auto query_span = parse_span(span);
  const auto table = std::visit(
      [&](const auto& sel_ptr) -> std::shared_ptr<arrow::Table> {
        assert(!!sel_ptr);
        return std::visit(
            [&]([[maybe_unused]] auto count) -> std::shared_ptr<arrow::Table> {
              using N = decltype(count);
              if (pixel_format == PixelFormat::BG2) {
                return make_bg2_arrow_df<N>(*sel_ptr, query_span);
              }
              assert(pixel_format == PixelFormat::COO);
              return make_coo_arrow_df<N>(*sel_ptr, query_span);
            },
            pixel_count);
      },
      selector);

  return nb::steal(arrow::py::wrap_table(table));
}

nb::object PixelSelector::to_pandas(std::string_view span) const {
  return to_arrow(span).attr("to_pandas")();
}

nb::object PixelSelector::to_df(std::string_view span) const { return to_pandas(span); }

template <typename N, typename PixelSelector>
[[nodiscard]] static nb::object make_csr_matrix(std::shared_ptr<const PixelSelector> sel,
                                                hictk::transformers::QuerySpan span) {
  if constexpr (std::is_same_v<N, long double>) {
    return make_csr_matrix<double>(std::move(sel), span);
  } else {
    return nb::cast(hictk::transformers::ToSparseMatrix(std::move(sel), N{}, span)());
  }
}

nb::object PixelSelector::to_csr(std::string_view span) const {
  const auto query_span = parse_span(span);

  return std::visit(
      [&](auto sel_ptr) -> nb::object {
        return std::visit(
            [&]([[maybe_unused]] auto count) -> nb::object {
              using N = decltype(count);
              return make_csr_matrix<N>(std::move(sel_ptr), query_span);
            },
            pixel_count);
      },
      selector);
}

nb::object PixelSelector::to_coo(std::string_view span) const {
  return to_csr(span).attr("tocoo")(false);
}

template <typename N, typename PixelSelector>
[[nodiscard]] static nb::object make_numpy_matrix(std::shared_ptr<const PixelSelector> sel,
                                                  hictk::transformers::QuerySpan span) {
  if constexpr (std::is_same_v<N, long double>) {
    return make_numpy_matrix<double>(std::move(sel), span);
  } else {
    return nb::cast(hictk::transformers::ToDenseMatrix(std::move(sel), N{}, span)());
  }
}

nb::object PixelSelector::to_numpy(std::string_view span) const {
  const auto query_span = parse_span(span);

  return std::visit(
      [&](auto sel_ptr) -> nb::object {
        return std::visit(
            [&]([[maybe_unused]] auto count) -> nb::object {
              using N = decltype(count);
              return make_numpy_matrix<N>(std::move(sel_ptr), query_span);
            },
            pixel_count);
      },
      selector);
}

template <typename N, typename PixelSelector>
[[nodiscard]] static nb::object compute_pixel_sum(const PixelSelector& sel) {
  return nb::cast(std::accumulate(
      sel.template begin<N>(), sel.template end<N>(), N{0},
      [](N accumulator, const hictk::ThinPixel<N>& tp) { return accumulator + tp.count; }));
}

nb::object PixelSelector::sum() const {
  const bool fp_sum = std::visit(
      [&]([[maybe_unused]] const auto count) { return std::is_floating_point_v<decltype(count)>; },
      pixel_count);

  const bool unsigned_sum = std::visit(
      [&]([[maybe_unused]] const auto count) { return std::is_unsigned_v<decltype(count)>; },
      pixel_count);

  return std::visit(
      [&](const auto& sel_ptr) -> nb::object {
        if (fp_sum) {
          return compute_pixel_sum<double>(*sel_ptr);
        }
        if (unsigned_sum) {
          return compute_pixel_sum<std::uint64_t>(*sel_ptr);
        }
        return compute_pixel_sum<std::int64_t>(*sel_ptr);
      },
      selector);
}

std::int64_t PixelSelector::nnz() const {
  return std::visit(
      [&](const auto& s) {
        using T = std::int_fast8_t;
        return std::distance(s->template begin<T>(), s->template end<T>());
      },
      selector);
}

hictk::transformers::QuerySpan PixelSelector::parse_span(std::string_view span) {
  if (span == "upper_triangle") {
    return hictk::transformers::QuerySpan::upper_triangle;
  }
  if (span == "lower_triangle") {
    return hictk::transformers::QuerySpan::lower_triangle;
  }
  if (span == "full") {
    return hictk::transformers::QuerySpan::full;
  }

  throw std::runtime_error(
      fmt::format(FMT_STRING("unrecognized query span \"{}\". Supported query spans are: "
                             "upper_triangle, lower_triangle, and full"),
                  span));
}

hictk::internal::NumericVariant PixelSelector::parse_count_type(std::string_view type) {
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

  sel.def(nb::init<std::shared_ptr<const hictk::cooler::PixelSelector>, std::string_view, bool>(),
          nb::arg("selector"), nb::arg("type"), nb::arg("join"));
  sel.def(nb::init<std::shared_ptr<const hictk::hic::PixelSelector>, std::string_view, bool>(),
          nb::arg("selector"), nb::arg("type"), nb::arg("join"));
  sel.def(nb::init<std::shared_ptr<const hictk::hic::PixelSelectorAll>, std::string_view, bool>(),
          nb::arg("selector"), nb::arg("type"), nb::arg("join"));

  sel.def("__repr__", &PixelSelector::repr);

  sel.def("coord1", &PixelSelector::get_coord1, "Get query coordinates for the first dimension.");
  sel.def("coord2", &PixelSelector::get_coord2, "Get query coordinates for the second dimension.");

  sel.def("__iter__", &PixelSelector::make_iterable, nb::keep_alive<0, 1>(),
          nb::sig("def __iter__(self) -> PixelIterator"),
          "Return an iterator over the selected pixels.");

  sel.def("to_arrow", &PixelSelector::to_arrow, nb::arg("query_span") = "upper_triangle",
          nb::sig("def to_arrow(self, query_span: str = \"upper_triangle\") -> pyarrow.Table"),
          "Retrieve interactions as a pandas DataFrame.");
  sel.def("to_pandas", &PixelSelector::to_pandas, nb::arg("query_span") = "upper_triangle",
          nb::sig("def to_pandas(self, query_span: str = \"upper_triangle\") -> pandas.DataFrame"),
          "Retrieve interactions as a pandas DataFrame.");
  sel.def("to_df", &PixelSelector::to_df, nb::arg("query_span") = "upper_triangle",
          nb::sig("def to_df(self, query_span: str = \"upper_triangle\") -> pandas.DataFrame"),
          "Alias to to_pandas().");
  sel.def("to_numpy", &PixelSelector::to_numpy, nb::arg("query_span") = "full",
          nb::sig("def to_numpy(self, query_span: str = \"full\") -> numpy.ndarray"),
          "Retrieve interactions as a numpy 2D matrix.");
  sel.def(
      "to_coo", &PixelSelector::to_coo, nb::arg("query_span") = "upper_triangle",
      nb::sig("def to_coo(self, query_span: str = \"upper_triangle\") -> scipy.sparse.coo_matrix"),
      "Retrieve interactions as a SciPy COO matrix.");
  sel.def(
      "to_csr", &PixelSelector::to_csr, nb::arg("query_span") = "upper_triangle",
      nb::sig("def to_csr(self, query_span: str = \"upper_triangle\") -> scipy.sparse.csr_matrix"),
      "Retrieve interactions as a SciPy CSR matrix.");

  sel.def("nnz", &PixelSelector::nnz,
          "Get the number of non-zero entries for the current pixel selection.");
  sel.def("sum", &PixelSelector::sum, nb::sig("def sum(self) -> int | float"),
          "Get the total number of interactions for the current pixel selection.");
}

}  // namespace hictkpy