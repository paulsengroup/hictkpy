// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <pybind11/pybind11.h>

#include <variant>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/pixel_selector.hpp"

namespace py = pybind11;

namespace hictkpy {

PixelSelector::PixelSelector(std::shared_ptr<const hictk::cooler::PixelSelector> sel_,
                             std::string_view type, bool join_)
    : selector(std::move(sel_)), join(join_) {
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
                             std::string_view type, bool join_)
    : selector(std::move(sel_)), join(join_) {
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
                             std::string_view type, bool join_)
    : selector(std::move(sel_)), join(join_) {
  if (type != "int" && type != "float") {
    throw std::runtime_error("type should be int or float");
  }

  if (type == "int") {
    pixel_count = std::int32_t{};
  } else {
    pixel_count = double{};
  }
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
  const auto c = coord1();
  return PixelCoordTuple{std::make_tuple(c.bin1.chrom().name(), c.bin1.start(), c.bin1.end(),
                                         c.bin2.chrom().name(), c.bin2.start(), c.bin2.end())};
}

auto PixelSelector::get_coord2() const -> PixelCoordTuple {
  const auto c = coord2();
  return PixelCoordTuple{std::make_tuple(c.bin1.chrom().name(), c.bin1.start(), c.bin1.end(),
                                         c.bin2.chrom().name(), c.bin2.start(), c.bin2.end())};
}

py::iterator PixelSelector::make_iterable() const {
  if (join) {
    return std::visit(
        [&](const auto& s) {
          if (int_pixels()) {
            using T = std::int32_t;
            auto jsel = hictk::transformers::JoinGenomicCoords(
                s->template begin<T>(), s->template end<T>(),
                std::make_shared<const hictk::BinTable>(bins()));
            return py::make_iterator(jsel.begin(), jsel.end());
          }
          using T = double;
          auto jsel = hictk::transformers::JoinGenomicCoords(
              s->template begin<T>(), s->template end<T>(),
              std::make_shared<const hictk::BinTable>(bins()));
          return py::make_iterator(jsel.begin(), jsel.end());
        },
        selector);
  }
  return std::visit(
      [&](const auto& s) {
        if (int_pixels()) {
          using T = std::int32_t;
          return py::make_iterator(s->template begin<T>(), s->template end<T>());
        }
        using T = double;
        return py::make_iterator(s->template begin<T>(), s->template end<T>());
      },
      selector);
}

py::object PixelSelector::to_df() const {
  return std::visit(
      [&](const auto& s) {
        if (int_pixels()) {
          using T = std::int32_t;
          return pixel_iterators_to_df(s->bins(), s->template begin<T>(), s->template end<T>(),
                                       join);
        } else {
          using T = double;
          return pixel_iterators_to_df(s->bins(), s->template begin<T>(), s->template end<T>(),
                                       join);
        }
      },
      selector);
}

py::object PixelSelector::to_coo() const {
  const auto bin_size = bins().bin_size();

  const auto span1 = coord1().bin2.end() - coord1().bin1.start();
  const auto span2 = coord2().bin2.end() - coord2().bin1.start();
  const auto num_rows = span1 == 0 ? bins().size() : (span1 + bin_size - 1) / bin_size;
  const auto num_cols = span2 == 0 ? bins().size() : (span2 + bin_size - 1) / bin_size;
  return std::visit(
      [&, num_rows, num_cols](const auto& s) {
        if (int_pixels()) {
          using T = std::int32_t;
          return pixel_iterators_to_coo(s->template begin<T>(), s->template end<T>(), num_rows,
                                        num_cols, coord1().bin1.id(), coord2().bin1.id());
        } else {
          using T = double;
          return pixel_iterators_to_coo(s->template begin<T>(), s->template end<T>(), num_rows,
                                        num_cols, coord1().bin1.id(), coord2().bin1.id());
        }
      },
      selector);
}

py::object PixelSelector::to_numpy() const {
  const auto bin_size = bins().bin_size();

  const auto span1 = coord1().bin2.end() - coord1().bin1.start();
  const auto span2 = coord2().bin2.end() - coord2().bin1.start();
  const auto num_rows = span1 == 0 ? bins().size() : (span1 + bin_size - 1) / bin_size;
  const auto num_cols = span2 == 0 ? bins().size() : (span2 + bin_size - 1) / bin_size;

  const auto mirror_matrix = coord1().bin1.chrom() == coord2().bin1.chrom();

  return std::visit(
      [&, num_rows, num_cols, mirror_matrix](const auto& s) {
        if (int_pixels()) {
          using T = std::int32_t;
          return pixel_iterators_to_numpy(s->template begin<T>(), s->template end<T>(), num_rows,
                                          num_cols, mirror_matrix, coord1().bin1.id(),
                                          coord2().bin1.id());
        } else {
          using T = double;
          return pixel_iterators_to_numpy(s->template begin<T>(), s->template end<T>(), num_rows,
                                          num_cols, mirror_matrix, coord1().bin1.id(),
                                          coord2().bin1.id());
        }
      },
      selector);
}

py::object PixelSelector::sum() const {
  return std::visit(
      [&](const auto& s) -> py::object {
        if (int_pixels()) {
          using T = std::int32_t;
          return py::cast(
              std::accumulate(s->template begin<T>(), s->template end<T>(), std::int64_t(0),
                              [](std::int64_t accumulator, const hictk::ThinPixel<T>& tp) {
                                return accumulator + tp.count;
                              }));
        } else {
          using T = double;
          return py::cast(std::accumulate(
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
