// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/pixel.hpp"

#include <fmt/compile.h>
#include <fmt/format.h>

#include <cstddef>
#include <cstdint>
#include <hictk/bin.hpp>
#include <hictk/pixel.hpp>
#include <string>
#include <string_view>
#include <variant>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/type.hpp"
#include "hictkpy/variant.hpp"

namespace hictkpy {

// NOLINTNEXTLINE(bugprone-exception-escape)
std::variant<std::int64_t, double> Pixel::count() const noexcept {
  return std::visit(
      [&](const auto n) -> std::variant<std::int64_t, double> {
        using N = remove_cvref_t<decltype(n)>;
        if constexpr (std::is_integral_v<N>) {
          return conditional_static_cast<std::int64_t>(n);
        } else {
          return conditional_static_cast<double>(n);
        }
      },
      _count);
}

const hictk::PixelCoordinates &Pixel::coords() const {
  if (HICTKPY_LIKELY(_coords.has_value())) {
    return *_coords;
  }

  throw nanobind::attribute_error(
      "Pixel does not have Bin with genomic coordinates associated with it. "
      "If you intend to access the genomic coordinates of Pixels, please make sure to call "
      "PixelSelector.fetch() with join=True.");
}

const hictk::Bin &Pixel::bin1() const { return coords().bin1; }
const hictk::Bin &Pixel::bin2() const { return coords().bin2; }

std::string_view Pixel::chrom1() const { return bin1().chrom().name(); }
std::int64_t Pixel::start1() const { return static_cast<std::int64_t>(bin1().start()); }
std::int64_t Pixel::end1() const { return static_cast<std::int64_t>(bin1().end()); }
std::string_view Pixel::chrom2() const { return bin2().chrom().name(); }
std::int64_t Pixel::start2() const { return static_cast<std::int64_t>(bin2().start()); }
std::int64_t Pixel::end2() const { return static_cast<std::int64_t>(bin2().end()); }

std::string Pixel::repr() const {
  return std::visit(
      [&](const auto &n) {
        if (_coords.has_value()) {
          return fmt::format(
              FMT_COMPILE("chrom1={}; start1={}; end1={}; chrom2={}; start2={}; end2={};"),
              _coords->bin1.chrom().name(), _coords->bin1.start(), _coords->bin1.end(),
              _coords->bin2.chrom().name(), _coords->bin2.start(), _coords->bin2.end(), n);
        }

        return fmt::format(FMT_COMPILE("bin1_id={}; bin2_id={}; count={};"), _bin1_id, _bin2_id, n);
      },
      _count);
}

std::string Pixel::str() const {
  return std::visit(
      [&](const auto &n) {
        if (_coords.has_value()) {
          return fmt::format(FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{}"), _coords->bin1.chrom().name(),
                             _coords->bin1.start(), _coords->bin1.end(),
                             _coords->bin2.chrom().name(), _coords->bin2.start(),
                             _coords->bin2.end(), n);
        }
        return fmt::format(FMT_COMPILE("{}\t{}\t{}"), _bin1_id, _bin2_id, n);
      },
      _count);
}

template <std::size_t I>
static void register_pixel_class_helper(nanobind::module_ &m, nanobind::class_<Pixel> &c) {
  using Var = NumericDtype;
  if constexpr (I < std::variant_size_v<Var>) {
    using N = std::variant_alternative_t<I, Var>;
    c.def(nanobind::init_implicit<const hictk::ThinPixel<N> &>(), "Private constructor.");
    c.def(nanobind::init_implicit<hictk::Pixel<N>>(), "Private constructor.");
    c.def(nanobind::init<hictk::Bin, hictk::Bin, N>(),
          "Construct a Pixel given a pair of Bins and the number of interactions.");
    c.def(nanobind::init<std::int64_t, std::int64_t, N>(),
          "Construct a Pixel given a pair of Bin identifiers and the number of interactions.");

    return register_pixel_class_helper<I + 1>(m, c);
  }
}

[[nodiscard]] static nanobind::class_<Pixel> register_pixel_class(nanobind::module_ &m) {
  auto pxl = nanobind::class_<Pixel>(m, "Pixel", "Class modeling a Pixel in COO or BG2 format.");
  register_pixel_class_helper<0>(m, pxl);

  return pxl;
}

void Pixel::bind(nanobind::module_ &m) {
  const auto *count_attr_sig = typing_union_required() ? "def count(self) -> Union[int, float]"
                                                       : "def count(self) -> int | float";

  register_pixel_class(m)
      .def_prop_ro("bin1_id", &Pixel::bin1_id, "Get the ID of bin1.")
      .def_prop_ro("bin2_id", &Pixel::bin2_id, "Get the ID of bin2.")
      .def_prop_ro("bin1", &Pixel::bin1, "Get bin1.", nanobind::rv_policy::copy)
      .def_prop_ro("bin2", &Pixel::bin2, "Get bin2.", nanobind::rv_policy::copy)
      .def_prop_ro("chrom1", &Pixel::chrom1, "Get the chromosome associated with bin1.")
      .def_prop_ro("start1", &Pixel::start1, "Get the start position associated with bin1.")
      .def_prop_ro("end1", &Pixel::end1, "Get the end position associated with bin1.")
      .def_prop_ro("chrom2", &Pixel::chrom2, "Get the chromosome associated with bin2.")
      .def_prop_ro("start2", &Pixel::start2, "Get the start position associated with bin2.")
      .def_prop_ro("end2", &Pixel::end2, "Get the end position associated with bin2.")
      .def_prop_ro("count", &Pixel::count, nanobind::sig(count_attr_sig),
                   "Get the number of interactions.")
      .def("__repr__", &Pixel::repr)
      .def("__str__", &Pixel::str);
}

}  // namespace hictkpy
