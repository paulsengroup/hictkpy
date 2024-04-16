// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>

#include <variant>

#include "hictk/cooler/cooler.hpp"
#include "hictk/fmt.hpp"
#include "hictk/hic.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/pixel_selector.hpp"

namespace nb = nanobind;

namespace hictkpy {

PixelSelector::PixelSelector(std::shared_ptr<const hictk::cooler::PixelSelector> sel_,
                             std::string_view type, bool join_, bool mirror_)
    : selector(std::move(sel_)), join(join_), mirror(mirror_) {
  if (type != "int" && type != "float") {
    throw std::runtime_error("type should be int or float");
  }

  if (type == "int") {
    pixel_count = std::int32_t{};
  } else {
    pixel_count = double{};
  }
}

PixelSelector::PixelSelector(std::shared_ptr<const hictk::hic::PixelSelector> sel_,
                             std::string_view type, bool join_, bool mirror_)
    : selector(std::move(sel_)), join(join_), mirror(mirror_) {
  if (type != "int" && type != "float") {
    throw std::runtime_error("type should be int or float");
  }

  if (type == "int") {
    pixel_count = std::int32_t{};
  } else {
    pixel_count = double{};
  }
}

PixelSelector::PixelSelector(std::shared_ptr<const hictk::hic::PixelSelectorAll> sel_,
                             std::string_view type, bool join_, bool mirror_)
    : selector(std::move(sel_)), join(join_), mirror(mirror_) {
  if (type != "int" && type != "float") {
    throw std::runtime_error("type should be int or float");
  }

  if (type == "int") {
    pixel_count = std::int32_t{};
  } else {
    pixel_count = double{};
  }
}

std::string PixelSelector::repr() const {
  if (!coord1()) {
    return "PixelSelector(ALL)";
  }

  return fmt::format(FMT_STRING("PixelSelector({}, {})"), coord1(), coord2());
}

constexpr bool PixelSelector::int_pixels() const noexcept {
  return std::holds_alternative<std::int32_t>(pixel_count);
}

constexpr bool PixelSelector::float_pixels() const noexcept { return !int_pixels(); }

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
  const auto c = mirror ? coord2() : coord1();
  return PixelCoordTuple{std::make_tuple(c.bin1.chrom().name(), c.bin1.start(), c.bin1.end(),
                                         c.bin2.chrom().name(), c.bin2.start(), c.bin2.end())};
}

auto PixelSelector::get_coord2() const -> PixelCoordTuple {
  const auto c = mirror ? coord1() : coord2();
  return PixelCoordTuple{std::make_tuple(c.bin1.chrom().name(), c.bin1.start(), c.bin1.end(),
                                         c.bin2.chrom().name(), c.bin2.start(), c.bin2.end())};
}

nb::iterator PixelSelector::make_iterable() const {
  if (mirror) {
    throw std::runtime_error(
        "iterating through the pixels for a query overlapping the lower triangle is not supported");
  }

  if (join) {
    return std::visit(
        [&](const auto& s) {
          if (int_pixels()) {
            using T = std::int32_t;
            auto jsel = hictk::transformers::JoinGenomicCoords(
                s->template begin<T>(), s->template end<T>(),
                std::make_shared<const hictk::BinTable>(bins()));
            return nb::make_iterator(nb::type<PixelSelector>(), "PixelIterator", jsel.begin(),
                                     jsel.end());
          }
          using T = double;
          auto jsel = hictk::transformers::JoinGenomicCoords(
              s->template begin<T>(), s->template end<T>(),
              std::make_shared<const hictk::BinTable>(bins()));
          return nb::make_iterator(nb::type<PixelSelector>(), "PixelIterator", jsel.begin(),
                                   jsel.end());
        },
        selector);
  }
  return std::visit(
      [&](const auto& s) {
        if (int_pixels()) {
          using T = std::int32_t;
          return nb::make_iterator(nb::type<PixelSelector>(), "PixelIterator",
                                   s->template begin<T>(), s->template end<T>());
        }
        using T = double;
        return nb::make_iterator(nb::type<PixelSelector>(), "PixelIterator", s->template begin<T>(),
                                 s->template end<T>());
      },
      selector);
}

nb::object PixelSelector::to_df() const {
  return std::visit(
      [&](const auto& s) {
        if (int_pixels()) {
          using T = std::int32_t;
          return pixel_iterators_to_df(s->bins(), s->template begin<T>(), s->template end<T>(),
                                       join, mirror);
        } else {
          using T = double;
          return pixel_iterators_to_df(s->bins(), s->template begin<T>(), s->template end<T>(),
                                       join, mirror);
        }
      },
      selector);
}

nb::object PixelSelector::to_coo() const {
  const auto bin_size = bins().resolution();

  const auto span1 = coord1().bin2.end() - coord1().bin1.start();
  const auto span2 = coord2().bin2.end() - coord2().bin1.start();
  const auto num_rows = span1 == 0 ? bins().size() : (span1 + bin_size - 1) / bin_size;
  const auto num_cols = span2 == 0 ? bins().size() : (span2 + bin_size - 1) / bin_size;
  return std::visit(
      [&, num_rows, num_cols](const auto& s) {
        if (int_pixels()) {
          using T = std::int32_t;
          return pixel_iterators_to_coo(s->template begin<T>(), s->template end<T>(), num_rows,
                                        num_cols, mirror, coord1().bin1.id(), coord2().bin1.id());
        } else {
          using T = double;
          return pixel_iterators_to_coo(s->template begin<T>(), s->template end<T>(), num_rows,
                                        num_cols, mirror, coord1().bin1.id(), coord2().bin1.id());
        }
      },
      selector);
}

nb::object PixelSelector::to_numpy() const {
  const auto bin_size = bins().resolution();

  const auto span1 = coord1().bin2.end() - coord1().bin1.start();
  const auto span2 = coord2().bin2.end() - coord2().bin1.start();
  const auto num_rows = span1 == 0 ? bins().size() : (span1 + bin_size - 1) / bin_size;
  const auto num_cols = span2 == 0 ? bins().size() : (span2 + bin_size - 1) / bin_size;

  const auto mirror_matrix = coord1().bin1.chrom() == coord2().bin1.chrom();

  const auto row_offset = static_cast<std::int64_t>(coord1().bin1.id());
  const auto col_offset = static_cast<std::int64_t>(coord2().bin1.id());

  return std::visit(
      [&, num_rows, num_cols, mirror_matrix, row_offset, col_offset](const auto& s) {
        if (int_pixels()) {
          using T = std::int32_t;
          return pixel_iterators_to_numpy(s->template begin<T>(), s->template end<T>(), num_rows,
                                          num_cols, mirror_matrix, row_offset, col_offset);
        } else {
          using T = double;
          return pixel_iterators_to_numpy(s->template begin<T>(), s->template end<T>(), num_rows,
                                          num_cols, mirror_matrix, row_offset, col_offset);
        }
      },
      selector);
}

nb::object PixelSelector::sum() const {
  return std::visit(
      [&](const auto& s) -> nb::object {
        if (int_pixels()) {
          using T = std::int32_t;
          return nb::cast(
              std::accumulate(s->template begin<T>(), s->template end<T>(), std::int64_t(0),
                              [](std::int64_t accumulator, const hictk::ThinPixel<T>& tp) {
                                return accumulator + tp.count;
                              }));
        } else {
          using T = double;
          return nb::cast(std::accumulate(
              s->template begin<T>(), s->template end<T>(), double(0),
              [](T accumulator, const hictk::ThinPixel<T>& tp) { return accumulator + tp.count; }));
        }
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

}  // namespace hictkpy
