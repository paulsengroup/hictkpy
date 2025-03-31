// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <hictk/pixel.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

template <typename N>
inline void declare_thin_pixel_class(nanobind::module_ &m, const std::string &suffix) {
  const auto type_name = std::string{"ThinPixel"} + suffix;
  nanobind::class_<hictk::ThinPixel<N>>(m, type_name.c_str(), "Pixel in COO format.")
      .def_prop_ro("bin1_id", [](const hictk::ThinPixel<N> &tp) { return tp.bin1_id; })
      .def_prop_ro("bin2_id", [](const hictk::ThinPixel<N> &tp) { return tp.bin2_id; })
      .def_prop_ro("count", [](const hictk::ThinPixel<N> &tp) { return tp.count; })
      .def("__repr__",
           [](const hictk::ThinPixel<N> &tp) {
             return fmt::format(FMT_COMPILE("bin1_id={}; bin2_id={}; count={};"), tp.bin1_id,
                                tp.bin2_id, tp.count);
           })
      .def("__str__", [](const hictk::ThinPixel<N> &tp) {
        return fmt::format(FMT_COMPILE("{}\t{}\t{}"), tp.bin1_id, tp.bin2_id, tp.count);
      });
}

template <typename N>
inline void declare_pixel_class(nanobind::module_ &m, const std::string &suffix) {
  const auto type_name = std::string{"Pixel"} + suffix;
  nanobind::class_<hictk::Pixel<N>>(m, type_name.c_str(), "Pixel in BG2 format.")
      .def_prop_ro("bin1_id", [](const hictk::Pixel<N> &p) { return p.coords.bin1.id(); })
      .def_prop_ro("bin2_id", [](const hictk::Pixel<N> &p) { return p.coords.bin2.id(); })
      .def_prop_ro("rel_bin1_id", [](const hictk::Pixel<N> &p) { return p.coords.bin1.rel_id(); })
      .def_prop_ro("rel_bin2_id", [](const hictk::Pixel<N> &p) { return p.coords.bin2.rel_id(); })
      .def_prop_ro("chrom1", [](const hictk::Pixel<N> &p) { return p.coords.bin1.chrom().name(); })
      .def_prop_ro("start1", [](const hictk::Pixel<N> &p) { return p.coords.bin1.start(); })
      .def_prop_ro("end1", [](const hictk::Pixel<N> &p) { return p.coords.bin1.end(); })
      .def_prop_ro("chrom2", [](const hictk::Pixel<N> &p) { return p.coords.bin2.chrom().name(); })
      .def_prop_ro("start2", [](const hictk::Pixel<N> &p) { return p.coords.bin2.start(); })
      .def_prop_ro("end2", [](const hictk::Pixel<N> &p) { return p.coords.bin2.end(); })
      .def_prop_ro("count", [](const hictk::Pixel<N> &p) { return p.count; })
      .def("__repr__",
           [](const hictk::Pixel<N> &p) {
             return fmt::format(
                 FMT_COMPILE("chrom1={}; start1={}; end1={}; chrom2={}; start2={}; end2={};"),
                 p.coords.bin1.chrom().name(), p.coords.bin1.start(), p.coords.bin1.end(),
                 p.coords.bin2.chrom().name(), p.coords.bin2.start(), p.coords.bin2.end(), p.count);
           })
      .def("__str__", [](const hictk::Pixel<N> &p) {
        return fmt::format(FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{}"), p.coords.bin1.chrom().name(),
                           p.coords.bin1.start(), p.coords.bin1.end(), p.coords.bin2.chrom().name(),
                           p.coords.bin2.start(), p.coords.bin2.end(), p.count);
      });
}

template <typename N>
inline std::vector<hictk::ThinPixel<N>> coo_df_to_thin_pixels(nanobind::object df, bool sort) {
  using BufferT1 = nanobind::ndarray<nanobind::numpy, nanobind::shape<-1>, std::uint64_t>;
  using BufferT2 = nanobind::ndarray<nanobind::numpy, nanobind::shape<-1>, N>;

  auto bin1_ids_np = nanobind::cast<BufferT1>(df.attr("__getitem__")("bin1_id").attr("to_numpy")());
  auto bin2_ids_np = nanobind::cast<BufferT1>(df.attr("__getitem__")("bin2_id").attr("to_numpy")());
  auto counts_np = nanobind::cast<BufferT2>(df.attr("__getitem__")("count").attr("to_numpy")());

  const auto bin1_ids = bin1_ids_np.view();
  const auto bin2_ids = bin2_ids_np.view();
  const auto counts = counts_np.view();

  std::vector<hictk::ThinPixel<N>> buffer(bin1_ids_np.size());
  for (std::size_t i = 0; i < bin1_ids_np.size(); ++i) {
    buffer[i] = hictk::ThinPixel<N>{bin1_ids(i), bin2_ids(i), counts(i)};
  }

  if (sort) {
    std::sort(buffer.begin(), buffer.end());
  }

  return buffer;
}

template <typename N>
inline std::vector<hictk::ThinPixel<N>> bg2_df_to_thin_pixels(const hictk::BinTable &bin_table,
                                                              nanobind::object df, bool sort) {
  using BufferT1 = nanobind::ndarray<nanobind::numpy, nanobind::shape<-1>, std::uint32_t>;
  using BufferT2 = nanobind::ndarray<nanobind::numpy, nanobind::shape<-1>, N>;

  auto chrom1 = nanobind::cast<nanobind::list>(df.attr("__getitem__")("chrom1").attr("tolist")());
  auto start1_np = nanobind::cast<BufferT1>(df.attr("__getitem__")("start1").attr("to_numpy")());
  auto end1_np = nanobind::cast<BufferT1>(df.attr("__getitem__")("end1").attr("to_numpy")());
  auto chrom2 = nanobind::cast<nanobind::list>(df.attr("__getitem__")("chrom2").attr("tolist")());
  auto start2_np = nanobind::cast<BufferT1>(df.attr("__getitem__")("start2").attr("to_numpy")());
  auto end2_np = nanobind::cast<BufferT1>(df.attr("__getitem__")("end2").attr("to_numpy")());
  auto counts_np = nanobind::cast<BufferT2>(df.attr("__getitem__")("count").attr("to_numpy")());

  const auto start1 = start1_np.view();
  const auto end1 = end1_np.view();
  const auto start2 = start2_np.view();
  const auto end2 = end2_np.view();
  const auto counts = counts_np.view();

  const auto &reference = bin_table.chromosomes();
  std::vector<hictk::ThinPixel<N>> buffer(start1_np.size());
  for (std::size_t i = 0; i < start1_np.size(); ++i) {
    if (end1(i) < start1(i) || end2(i) < start2(i)) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Found an invalid pixel {} {} {} {} {} {} {}: bin end position cannot be "
                     "smaller than its start"),
          nanobind::cast<nanobind::str>(chrom1[i]).c_str(), start1(i), end1(i),
          nanobind::cast<nanobind::str>(chrom2[i]).c_str(), start2(i), end2(i), counts(i)));
    }

    auto bin1 =
        bin_table.at(reference.at(nanobind::cast<nanobind::str>(chrom1[i]).c_str()), start1(i));
    auto bin2 =
        bin_table.at(reference.at(nanobind::cast<nanobind::str>(chrom2[i]).c_str()), start2(i));

    if (bin_table.type() == hictk::BinTable::Type::fixed &&
        (end1(i) - start1(i) > bin_table.resolution() ||
         end2(i) - start2(i) > bin_table.resolution())) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Found an invalid pixel {} {} {} {} {} {} {}: pixel spans a "
                     "distance greater than the bin size"),
          nanobind::cast<nanobind::str>(chrom1[i]).c_str(), start1(i), end1(i),
          nanobind::cast<nanobind::str>(chrom2[i]).c_str(), start2(i), end2(i), counts(i)));
    }

    buffer[i] =
        hictk::Pixel<N>{hictk::PixelCoordinates{std::move(bin1), std::move(bin2)}, counts(i)}
            .to_thin();
  }

  if (sort) {
    std::sort(buffer.begin(), buffer.end());
  }

  return buffer;
}

}  // namespace hictkpy
