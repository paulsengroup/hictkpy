// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef _WIN32
// Workaround bug several symbol redefinition errors due to something including <winsock.h>
#include <winsock2.h>
#endif

#include <fmt/format.h>

// clang-format off
#include "hictkpy/suppress_warnings.hpp"
HICTKPY_DISABLE_WARNING_PUSH
HICTKPY_DISABLE_WARNING_OLD_STYLE_CAST
HICTKPY_DISABLE_WARNING_PEDANTIC
HICTKPY_DISABLE_WARNING_SHADOW
HICTKPY_DISABLE_WARNING_SIGN_CONVERSION
HICTKPY_DISABLE_WARNING_USELESS_CAST
#include <arrow/python/arrow_to_pandas.h>

#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
HICTKPY_DISABLE_WARNING_POP
// clang-format on

#include <variant>

#include "hictk/cooler/cooler.hpp"
#include "hictk/fmt.hpp"
#include "hictk/hic.hpp"
#include "hictk/transformers/to_dataframe.hpp"
#include "hictk/transformers/to_dense_matrix.hpp"
#include "hictk/transformers/to_sparse_matrix.hpp"
#include "hictkpy/common.hpp"
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

hictk::PixelCoordinates PixelSelector::coord1() const noexcept {
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

hictk::PixelCoordinates PixelSelector::coord2() const noexcept {
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

const hictk::BinTable& PixelSelector::bins() const noexcept {
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

nb::object PixelSelector::make_iterable() const {
  return std::visit(
      [&](const auto& sel_ptr) -> nb::object {
        assert(!!sel_ptr);
        return std::visit(
            [&]([[maybe_unused]] auto count) -> nb::object {
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
[[nodiscard]] std::shared_ptr<arrow::Table> make_bg2_arrow_df(const PixelSelector& sel,
                                                              hictk::transformers::QuerySpan span) {
  if constexpr (std::is_same_v<N, long double>) {
    return make_bg2_arrow_df<double>(sel, span);
  } else {
    return hictk::transformers::ToDataFrame(sel.template begin<N>(), sel.template end<N>(),
                                            hictk::transformers::DataFrameFormat::BG2,
                                            sel.bins_ptr(), span)();
  }
}

template <typename N, typename PixelSelector>
[[nodiscard]] std::shared_ptr<arrow::Table> make_coo_arrow_df(const PixelSelector& sel,
                                                              hictk::transformers::QuerySpan span) {
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
[[nodiscard]] nb::object make_csr_matrix(std::shared_ptr<const PixelSelector> sel,
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
[[nodiscard]] nb::object make_numpy_matrix(std::shared_ptr<const PixelSelector> sel,
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
  static_assert(sizeof(unsigned) == 4);
  static_assert(sizeof(int) == 4);
  static_assert(sizeof(float) == 4);
  static_assert(sizeof(double) == 8);

  if (type == "uint8") {
    return {std::uint8_t{}};
  }
  if (type == "uint16") {
    return {std::uint16_t{}};
  }
  if (type == "uint32" || type == "uint") {
    return {std::uint32_t{}};
  }
  if (type == "uint64") {
    return {std::uint64_t{}};
  }
  if (type == "int8") {
    return {std::int8_t{}};
  }
  if (type == "int16") {
    return {std::int16_t{}};
  }
  if (type == "int32" || type == "int") {
    return {std::int32_t{}};
  }
  if (type == "int64") {
    return {std::int64_t{}};
  }
  if (type == "float32") {
    return {float{}};
  }
  if (type == "float64" || type == "float" || type == "double") {
    return {double{}};
  }

  throw std::runtime_error(fmt::format(
      FMT_STRING(
          "unable to map \"{}\" to a valid pixel count type. Valid types are: uint, int, float, "
          "double, uint8, uint16, uint32, uint64, int8, int16, int32, int64, float32, and float64"),
      type));
}

std::string_view PixelSelector::count_type_to_str(const PixelVar& var) noexcept {
  static_assert(sizeof(float) == 4);
  static_assert(sizeof(double) == 8);
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

  HICTK_UNREACHABLE_CODE;
}

}  // namespace hictkpy
