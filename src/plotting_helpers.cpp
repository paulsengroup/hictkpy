// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/plotting_helpers.hpp"

#include <cassert>
#include <cstdint>
#include <hictk/balancing/methods.hpp>
#include <hictk/balancing/weights.hpp>
#include <hictk/bin_table.hpp>
#include <hictk/file.hpp>
#include <hictk/pixel.hpp>
#include <hictk/transformers/coarsen.hpp>
#include <hictk/transformers/common.hpp>
#include <hictk/transformers/to_dense_matrix.hpp>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <variant>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"

namespace hictkpy {

namespace nb = nanobind;

[[nodiscard]] static constexpr std::uint32_t compute_target_resolution(
    std::uint32_t src, std::uint32_t target) noexcept {
  assert(src < target);
  assert(src != 0);

  const auto cfx = target / src;
  if (cfx < 2) {
    return src;
  }

  return src * cfx;
}

template <typename PixelIt>
class MockSelector {
  PixelIt _first{};
  PixelIt _last{};

  hictk::PixelCoordinates _coord1{};
  hictk::PixelCoordinates _coord2{};
  std::shared_ptr<const hictk::BinTable> _bins{};

 public:
  MockSelector() = delete;
  MockSelector(PixelIt first, PixelIt last, hictk::PixelCoordinates coord1_,
               hictk::PixelCoordinates coord2_, std::shared_ptr<const hictk::BinTable> bins_)
      : _first(std::move(first)),
        _last(std::move(last)),
        _coord1(std::move(coord1_)),
        _coord2(std::move(coord2_)),
        _bins(std::move(bins_)) {}

  template <typename M>
  [[nodiscard]] PixelIt begin() const {
    return _first;
  }
  template <typename M>
  [[nodiscard]] PixelIt end() const {
    return _last;
  }

  [[nodiscard]] const hictk::PixelCoordinates& coord1() const noexcept { return _coord1; }
  [[nodiscard]] const hictk::PixelCoordinates& coord2() const noexcept { return _coord2; }
  [[nodiscard]] const hictk::BinTable& bins() const noexcept { return *_bins; }
  [[nodiscard]] std::shared_ptr<const hictk::BinTable> bins_ptr() const noexcept { return _bins; }
  [[noreturn]] MockSelector fetch([[maybe_unused]] const hictk::PixelCoordinates& coord1_,
                                  [[maybe_unused]] const hictk::PixelCoordinates& coord2_) const {
    throw std::logic_error("MockSelector::fetch(): not implemented");
  }
  [[nodiscard]] const hictk::balancing::Weights& weights() const {
    static const hictk::balancing::Weights w{};
    return w;
  }
};

[[nodiscard]] static hictk::PixelSelector fetch(const hictk::File& f,
                                                const std::optional<std::string_view>& range1,
                                                const std::optional<std::string_view>& range2,
                                                std::string_view normalization) {
  if (range1.has_value() && range2.has_value()) {
    return f.fetch(*range1, *range2, hictk::balancing::Method{normalization});
  }
  if (range1.has_value()) {
    return f.fetch(*range1, hictk::balancing::Method{normalization});
  }
  return f.fetch(hictk::balancing::Method{normalization});
}

nanobind::object coarsen(const hictk::File& f, std::uint32_t target_resolution,
                         std::optional<std::string_view> range1,
                         std::optional<std::string_view> range2,
                         std::optional<std::string_view> normalization, bool floating_point,
                         bool exact) {
  // TODO ensure query fully overlaps with the new resolution (i.e. the last bin is either the last
  // bin on a given chromosome, or a multiple of the target resolution)
  // TODO we probably need to fix this dynamically
  assert(f.resolution() != 0);
  assert(target_resolution >= f.resolution());
  // TODO problem, need to do something about slice_weights(*_sel)

  const auto target_res_is_multiple = target_resolution % f.resolution() == 0;

  if (exact && !target_res_is_multiple) {
    throw std::invalid_argument("target_resolution must be a multiple of file resolution");
  }

  if (!target_res_is_multiple) {
    return coarsen(f, compute_target_resolution(f.resolution(), target_resolution),
                   std::move(range1), std::move(range2), std::move(normalization), floating_point,
                   true);
  }

  if (normalization.has_value()) {
    floating_point = true;
  } else {
    normalization = "NONE";
  }

  if (range1.has_value() && !range2.has_value()) {
    range2 = range1;
  }

  const auto sel = fetch(f, range1, range2, *normalization);

  std::variant<std::int64_t, double> pixel_type{std::int64_t{}};
  if (floating_point) {
    pixel_type = double{};
  }

  return std::visit(
      [&](const auto& sel_) {
        return std::visit(
            [&](const auto& count_type) -> nb::object {
              using N = remove_cvref_t<decltype(count_type)>;
              const hictk::transformers::CoarsenPixels zoomified_sel(
                  sel_.template begin<N>(), sel_.template end<N>(), f.bins_ptr(),
                  target_resolution / f.resolution());

              MockSelector mock_sel{zoomified_sel.begin(), zoomified_sel.end(), sel.coord1(),
                                    sel.coord2(), f.bins_ptr()};

              return nb::cast(hictk::transformers::ToDenseMatrix(std::move(mock_sel), N{})());
            },
            pixel_type);
      },
      sel.get());
}

}  // namespace hictkpy
