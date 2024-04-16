// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>

#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/numeric_variant.hpp"
#include "hictk/pixel.hpp"
#include "hictk/suppress_warnings.hpp"

namespace hictkpy {

namespace nb = nanobind;
using namespace nb::literals;

template <typename T>
[[nodiscard]] constexpr std::string_view map_type_to_dtype() {
  if constexpr (std::is_unsigned_v<T>) {
    switch (sizeof(T)) {
      case 1:
        return "uint8";
      case 2:
        return "uint16";
      case 4:
        return "uint32";
      case 8:
        return "uint64";
    }
  }
  if constexpr (std::is_integral_v<T>) {
    switch (sizeof(T)) {
      case 1:
        return "int8";
      case 2:
        return "int16";
      case 4:
        return "int32";
      case 8:
        return "int64";
    }
  }
  switch (sizeof(T)) {
    case 2:
      return "float16";
    case 4:
      return "float32";
    case 8:
      return "float64";
  }

  throw std::runtime_error("Unable to infer numpy dtype.");
}

[[nodiscard]] inline hictk::internal::NumericVariant map_dtype_to_type(std::string_view dtype) {
  if (dtype == "uint8") {
    return std::uint8_t{};
  }
  if (dtype == "uint16") {
    return std::uint16_t{};
  }
  if (dtype == "uint32") {
    return std::uint32_t{};
  }
  if (dtype == "uint64") {
    return std::uint64_t{};
  }

  if (dtype == "int8") {
    return std::int8_t{};
  }
  if (dtype == "int16") {
    return std::int16_t{};
  }
  if (dtype == "int32") {
    return std::int32_t{};
  }
  if (dtype == "int64") {
    return std::int64_t{};
  }

  if (dtype == "float16") {
    return float{};
  }
  if (dtype == "float32") {
    return float{};
  }
  if (dtype == "float64") {
    return double{};
  }

  throw std::runtime_error("Unable to map dtype " + std::string{dtype} + " to a C++ type.");
}

template <typename T>
struct Dynamic1DA {
 private:
  using BufferT = nb::ndarray<nb::numpy, nb::shape<nb::any>, T>;
  using VectorT = decltype(std::declval<BufferT>().view());
  nb::object _np_array{};
  BufferT _buff{};
  VectorT _vector{};
  nb::object _dtype{};

  std::int64_t _size{};
  std::int64_t _capacity{};

 public:
  inline explicit Dynamic1DA(std::size_t size_ = 1000)
      : _capacity(static_cast<std::int64_t>(size_)) {
    const auto np = nb::module_::import_("numpy");
    _dtype = np.attr("dtype")(map_type_to_dtype<T>());
    _np_array = np.attr("empty")(size_, "dtype"_a = _dtype);
    _buff = nb::cast<BufferT>(_np_array);
    _vector = _buff.view();
  }

  inline void push_back(T x) {
    if (_capacity == _size) {
      grow();
    }
    _vector(_size++) = x;
  }

  inline void emplace_back(T &&x) {
    if (_capacity == _size) {
      grow();
    }
    _vector(_size++) = std::move(x);
  }
  inline void resize(std::int64_t new_size) {
    if (_capacity == new_size) {
      return;
    }
    auto np = nb::module_::import_("numpy");

    auto new_array = np.attr("empty")(new_size, "dtype"_a = _dtype);
    auto new_buff = nb::cast<BufferT>(new_array);
    auto new_vector = new_buff.view();

    _capacity = new_size;
    _size = std::min(_capacity, _size);
    std::copy(_vector.data(), _vector.data() + static_cast<std::size_t>(_size), new_vector.data());

    std::swap(new_array, _np_array);
    std::swap(new_buff, _buff);
    std::swap(new_vector, _vector);
  }
  inline void grow() { resize(_buff.size() * 2); }
  inline void shrink_to_fit() { resize(_size); }
  [[nodiscard]] auto operator()() -> BufferT {
    shrink_to_fit();
    return _buff;
  }
};

template <typename File>
inline nb::dict get_chromosomes_from_file(const File &f, bool include_all = false) {
  nb::dict py_chroms{};  // NOLINT
  for (const auto &chrom : f.chromosomes()) {
    if (!include_all && chrom.is_all()) {
      continue;
    }
    const std::string name{chrom.name()};
    py_chroms[name.c_str()] = chrom.size();
  }

  return py_chroms;
}

template <typename File>
inline nb::object get_bins_from_file(const File &f) {
  auto pd = nb::module_::import_("pandas");

  Dynamic1DA<std::int32_t> chrom_ids{};
  Dynamic1DA<std::int32_t> starts{};
  Dynamic1DA<std::int32_t> ends{};
  for (const auto &bin : f.bins()) {
    chrom_ids.push_back(static_cast<std::int32_t>(bin.chrom().id()));
    starts.push_back(static_cast<std::int32_t>(bin.start()));
    ends.push_back(static_cast<std::int32_t>(bin.end()));
  }

  std::vector<std::string> chrom_names{};
  std::transform(f.chromosomes().begin(), f.chromosomes().end(), std::back_inserter(chrom_names),
                 [&](const hictk::Chromosome &chrom) { return std::string{chrom.name()}; });

  nb::dict py_bins_dict{};  // NOLINT

  py_bins_dict["chrom"] =
      pd.attr("Categorical")
          .attr("from_codes")(chrom_ids(), "categories"_a = chrom_names, "validate"_a = false);
  py_bins_dict["start"] = pd.attr("Series")(starts(), "copy"_a = false);
  py_bins_dict["end"] = pd.attr("Series")(ends(), "copy"_a = false);

  auto df = pd.attr("DataFrame")(py_bins_dict, "copy"_a = false);
  return df;
}

template <typename PixelIt>
inline nb::object pixel_iterators_to_coo(PixelIt first_pixel, PixelIt last_pixel,
                                         std::size_t num_rows, std::size_t num_cols, bool mirror,
                                         std::size_t row_offset = 0, std::size_t col_offset = 0) {
  using N = decltype(first_pixel->count);
  auto ss = nb::module_::import_("scipy.sparse");

  Dynamic1DA<std::int64_t> bin1_ids{};
  Dynamic1DA<std::int64_t> bin2_ids{};
  Dynamic1DA<N> counts{};

  std::for_each(first_pixel, last_pixel, [&](const hictk::ThinPixel<N> &tp) {
    bin1_ids.push_back(static_cast<std::int64_t>(tp.bin1_id - row_offset));
    bin2_ids.push_back(static_cast<std::int64_t>(tp.bin2_id - col_offset));
    counts.push_back(tp.count);
  });

  if (mirror) {
    std::swap(num_rows, num_cols);
    std::swap(bin1_ids, bin2_ids);
  }

  // See
  // https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html#scipy.sparse.coo_matrix
  // Building a sparse COO from an array triplet is much faster than converting
  // an Eigen matrix

  nb::list shape{};
  shape.append(num_rows);
  shape.append(num_cols);

  nb::list coords{};
  coords.append(bin1_ids());
  coords.append(bin2_ids());

  nb::list data{};
  data.append(counts());
  data.append(nb::tuple(coords));

  auto m = ss.attr("coo_matrix")(nb::tuple(data), "shape"_a = shape);
  return m;
}

template <typename PixelIt>
inline nb::object pixel_iterators_to_coo_df(PixelIt first_pixel, PixelIt last_pixel, bool mirror) {
  using N = decltype(first_pixel->count);

  auto pd = nb::module_::import_("pandas");

  Dynamic1DA<std::int64_t> bin1_ids{};
  Dynamic1DA<std::int64_t> bin2_ids{};
  Dynamic1DA<N> counts{};

  std::for_each(first_pixel, last_pixel, [&](const hictk::ThinPixel<N> &tp) {
    bin1_ids.push_back(static_cast<std::int64_t>(tp.bin1_id));
    bin2_ids.push_back(static_cast<std::int64_t>(tp.bin2_id));
    counts.push_back(tp.count);
  });

  if (mirror) {
    std::swap(bin1_ids, bin2_ids);
  }

  nb::dict py_pixels_dict{};  // NOLINT

  py_pixels_dict["bin1_id"] = pd.attr("Series")(bin1_ids(), "copy"_a = false);
  py_pixels_dict["bin2_id"] = pd.attr("Series")(bin2_ids(), "copy"_a = false);
  py_pixels_dict["count"] = pd.attr("Series")(counts(), "copy"_a = false);

  auto df = pd.attr("DataFrame")(py_pixels_dict, "copy"_a = false);

  if (mirror) {
    df.attr("sort_values")(std::vector<std::string>{"bin1_id", "bin2_id"}, "inplace"_a = true);
    df.attr("reset_index")();
  }

  return df;
}

template <typename PixelIt>
inline nb::object pixel_iterators_to_numpy(PixelIt first_pixel, PixelIt last_pixel,
                                           std::int64_t num_rows, std::int64_t num_cols,
                                           bool mirror_below_diagonal = true,
                                           std::size_t row_offset = 0, std::size_t col_offset = 0) {
  using N = decltype(first_pixel->count);

  const auto np = nb::module_::import_("numpy");
  const auto dtype = np.attr("dtype")(map_type_to_dtype<N>());
  auto buffer = np.attr("zeros")(std::vector<std::int64_t>{num_rows, num_cols}, "dtype"_a = dtype);

  using Shape = nb::shape<nb::any, nb::any>;
  using MatrixT = nb::ndarray<nb::numpy, N, Shape>;

  auto matrix = nb::cast<MatrixT>(buffer);
  auto m = matrix.view();

  std::for_each(first_pixel, last_pixel, [&](const hictk::ThinPixel<N> &tp) {
    const auto i1 = static_cast<std::int64_t>(tp.bin1_id - row_offset);
    const auto i2 = static_cast<std::int64_t>(tp.bin2_id - col_offset);
    m(i1, i2) = tp.count;

    if (mirror_below_diagonal) {
      const auto delta = i2 - i1;
      if (delta >= 0 && delta < num_rows && i1 < num_cols && i2 < num_rows) {
        m(i2, i1) = tp.count;
      } else if ((delta < 0 || delta > num_cols) && i1 < num_cols && i2 < num_rows) {
        const auto i3 = static_cast<std::int64_t>(tp.bin2_id - row_offset);
        const auto i4 = static_cast<std::int64_t>(tp.bin1_id - col_offset);

        if (i3 >= 0 && i3 < num_rows && i4 >= 0 && i4 < num_cols) {
          m(i3, i4) = tp.count;
        }
      }
    }
  });
  return nb::cast(matrix);
}

template <typename PixelIt>
inline nb::object pixel_iterators_to_bg2(const hictk::BinTable &bins, PixelIt first_pixel,
                                         PixelIt last_pixel, bool mirror) {
  using N = decltype(first_pixel->count);

  auto pd = nb::module_::import_("pandas");

  nb::ndarray<nb::numpy, int> test;

  Dynamic1DA<uint32_t> chrom1_ids{};
  Dynamic1DA<std::int32_t> starts1{};
  Dynamic1DA<std::int32_t> ends1{};
  Dynamic1DA<uint32_t> chrom2_ids{};
  Dynamic1DA<std::int32_t> starts2{};
  Dynamic1DA<std::int32_t> ends2{};
  Dynamic1DA<N> counts{};

  std::for_each(first_pixel, last_pixel, [&](const hictk::ThinPixel<N> &tp) {
    const hictk::Pixel<N> p{bins, tp};

    chrom1_ids.push_back(static_cast<std::int32_t>(p.coords.bin1.chrom().id()));
    starts1.push_back(static_cast<std::int32_t>(p.coords.bin1.start()));
    ends1.push_back(static_cast<std::int32_t>(p.coords.bin1.end()));

    chrom2_ids.push_back(static_cast<std::int32_t>(p.coords.bin2.chrom().id()));
    starts2.push_back(static_cast<std::int32_t>(p.coords.bin2.start()));
    ends2.push_back(static_cast<std::int32_t>(p.coords.bin2.end()));

    counts.push_back(p.count);
  });

  if (mirror) {
    std::swap(chrom1_ids, chrom2_ids);
    std::swap(starts1, starts2);
    std::swap(ends1, ends2);
  }

  std::vector<std::string> chrom_names{};
  std::transform(bins.chromosomes().begin(), bins.chromosomes().end(),
                 std::back_inserter(chrom_names),
                 [&](const hictk::Chromosome &chrom) { return std::string{chrom.name()}; });

  nb::dict py_pixels_dict{};  // NOLINT

  py_pixels_dict["chrom1"] =
      pd.attr("Categorical")
          .attr("from_codes")(chrom1_ids(), "categories"_a = chrom_names, "validate"_a = false);
  py_pixels_dict["start1"] = pd.attr("Series")(starts1(), "copy"_a = false);
  py_pixels_dict["end1"] = pd.attr("Series")(ends1(), "copy"_a = false);

  py_pixels_dict["chrom2"] =
      pd.attr("Categorical")
          .attr("from_codes")(chrom2_ids(), "categories"_a = chrom_names, "validate"_a = false);
  py_pixels_dict["start2"] = pd.attr("Series")(starts2(), "copy"_a = false);
  py_pixels_dict["end2"] = pd.attr("Series")(ends2(), "copy"_a = false);

  py_pixels_dict["count"] = pd.attr("Series")(counts(), "copy"_a = false);

  auto df = pd.attr("DataFrame")(py_pixels_dict, "copy"_a = false);
  if (mirror) {
    df.attr("sort_values")(std::vector<std::string>{"chrom1", "start1", "chrom2", "start2"},
                           "inplace"_a = true);
    df.attr("reset_index")();
  }

  return df;
}

template <typename PixelIt>
static nb::object pixel_iterators_to_df(const hictk::BinTable &bins, PixelIt first_pixel,
                                        PixelIt last_pixel, bool join, bool mirror) {
  if (join) {
    return pixel_iterators_to_bg2(bins, first_pixel, last_pixel, mirror);
  }
  return pixel_iterators_to_coo_df(first_pixel, last_pixel, mirror);
}

template <typename File>
inline nb::object file_fetch_all(File &f, std::string_view normalization,
                                 std::string_view count_type, bool join) {
  if (count_type != "int" && count_type != "float") {
    throw std::runtime_error("invalid count type. Allowed types: int, float.");
  }

  if (normalization != "NONE") {
    count_type = "float";
  }

  auto sel = f.fetch(hictk::balancing::Method{normalization});
  if (count_type == "int") {
    return pixel_iterators_to_df(f.bins(), sel.template begin<std::int32_t>(),
                                 sel.template end<std::int32_t>(), join);
  }
  return pixel_iterators_to_df(f.bins(), sel.template begin<double>(), sel.template end<double>(),
                               join);
}

template <typename File>
inline nb::object file_fetch(const File &f, std::string_view range1, std::string_view range2,
                             std::string_view normalization, std::string_view count_type, bool join,
                             std::string_view query_type) {
  if (range1.empty()) {
    return file_fetch_all(f, normalization, count_type, join);
  }
  if (normalization != "NONE") {
    count_type = "float";
  }

  const auto qt =
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED;

  auto sel = range2.empty() || range1 == range2
                 ? f.fetch(range1, hictk::balancing::Method(normalization), qt)
                 : f.fetch(range1, range2, hictk::balancing::Method(normalization), qt);

  if (count_type == "int") {
    return pixel_iterators_to_df(f.bins(), sel.template begin<std::int32_t>(),
                                 sel.template end<std::int32_t>(), join);
  }
  return pixel_iterators_to_df(f.bins(), sel.template begin<double>(), sel.template end<double>(),
                               join);
}

template <typename File>
inline nb::object file_fetch_all_sparse(File &f, std::string_view normalization,
                                        std::string_view count_type) {
  if (count_type != "int" && count_type != "float") {
    throw std::runtime_error("invalid count type. Allowed types: int, float.");
  }

  if (normalization != "NONE") {
    count_type = "float";
  }

  auto sel = f.fetch(hictk::balancing::Method{normalization});
  if (count_type == "int") {
    return pixel_iterators_to_coo(sel.template begin<std::int32_t>(),
                                  sel.template end<std::int32_t>(), f.bins().size(),
                                  f.bins().size());
  }
  return pixel_iterators_to_coo(sel.template begin<double>(), sel.template end<double>(),
                                f.bins().size(), f.bins().size());
}

template <typename File>
inline nb::object file_fetch_sparse(const File &f, std::string_view range1, std::string_view range2,
                                    std::string_view normalization, std::string_view count_type,
                                    std::string_view query_type) {
  if (range1.empty()) {
    return file_fetch_all_sparse(f, normalization, count_type);
  }
  if (normalization != "NONE") {
    count_type = "float";
  }

  const auto qt =
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED;

  const auto gi1 = hictk::GenomicInterval::parse(f.chromosomes(), std::string{range1}, qt);
  const auto gi2 = range2.empty()
                       ? gi1
                       : hictk::GenomicInterval::parse(f.chromosomes(), std::string{range2}, qt);

  const auto bin_size = f.bin_size();

  const auto num_rows = (gi1.size() + bin_size - 1) / bin_size;
  const auto num_cols = (gi2.size() + bin_size - 1) / bin_size;

  const auto bin1 = f.bins().at(gi1.chrom(), gi1.start());
  const auto bin2 = f.bins().at(gi2.chrom(), gi2.start());

  auto sel = range2.empty() || range1 == range2
                 ? f.fetch(range1, hictk::balancing::Method(normalization), qt)
                 : f.fetch(range1, range2, hictk::balancing::Method(normalization), qt);

  if (count_type == "int") {
    return pixel_iterators_to_coo(sel.template begin<std::int32_t>(),
                                  sel.template end<std::int32_t>(), num_rows, num_cols, bin1.id(),
                                  bin2.id());
  }
  return pixel_iterators_to_coo(sel.template begin<double>(), sel.template end<double>(), num_rows,
                                num_cols, bin1.id(), bin2.id());
}

template <typename File>
inline nb::object file_fetch_all_dense(File &f, std::string_view normalization,
                                       std::string_view count_type) {
  if (count_type != "int" && count_type != "float") {
    throw std::runtime_error("invalid count type. Allowed types: int, float.");
  }

  if (normalization != "NONE") {
    count_type = "float";
  }
  auto sel = f.fetch(hictk::balancing::Method{normalization});
  if (count_type == "int") {
    return pixel_iterators_to_numpy(sel.template begin<std::int32_t>(),
                                    sel.template end<std::int32_t>(), f.bins().size(),
                                    f.bins().size());
  }
  return pixel_iterators_to_numpy(sel.template begin<double>(), sel.template end<double>(),
                                  f.bins().size(), f.bins().size());
}

template <typename File>
inline nb::object file_fetch_dense(const File &f, std::string_view range1, std::string_view range2,
                                   std::string_view normalization, std::string_view count_type,
                                   std::string_view query_type) {
  if (range1.empty()) {
    return file_fetch_all_dense(f, normalization, count_type);
  }
  if (normalization != "NONE") {
    count_type = "float";
  }

  const auto qt =
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED;

  const auto gi1 = hictk::GenomicInterval::parse(f.chromosomes(), std::string{range1}, qt);
  const auto gi2 = range2.empty()
                       ? gi1
                       : hictk::GenomicInterval::parse(f.chromosomes(), std::string{range2}, qt);

  const auto bin_size = f.bin_size();

  const auto num_rows = (gi1.size() + bin_size - 1) / bin_size;
  const auto num_cols = (gi2.size() + bin_size - 1) / bin_size;

  const auto bin1 = f.bins().at(gi1.chrom(), gi1.start());
  const auto bin2 = f.bins().at(gi2.chrom(), gi2.start());

  const auto mirror_matrix = gi1.chrom() == gi2.chrom();

  auto sel = range2.empty() || range1 == range2
                 ? f.fetch(range1, hictk::balancing::Method(normalization), qt)
                 : f.fetch(range1, range2, hictk::balancing::Method(normalization), qt);

  if (count_type == "int") {
    return pixel_iterators_to_numpy(sel.template begin<std::int32_t>(),
                                    sel.template end<std::int32_t>(), num_rows, num_cols,
                                    mirror_matrix, bin1.id(), bin2.id());
  }
  return pixel_iterators_to_numpy(sel.template begin<double>(), sel.template end<double>(),
                                  num_rows, num_cols, mirror_matrix, bin1.id(), bin2.id());
}

template <typename File>
inline nb::object file_fetch_sum_all(const File &f, std::string_view normalization,
                                     std::string_view count_type) {
  if (normalization != "NONE") {
    count_type = "float";
  }

  auto sel = f.fetch(hictk::balancing::Method{normalization});
  if (count_type == "int") {
    return nb::cast(std::accumulate(
        sel.template begin<std::int32_t>(), sel.template end<std::int32_t>(), std::int64_t(0),
        [](const auto accumulator, const hictk::ThinPixel<std::int32_t> &p) {
          return accumulator + p.count;
        }));
  }
  return nb::cast(std::accumulate(sel.template begin<double>(), sel.template end<double>(), 0.0,
                                  [](const auto accumulator, const hictk::ThinPixel<double> &p) {
                                    return accumulator + p.count;
                                  }));
}

template <typename File>
inline nb::object file_fetch_sum(const File &f, std::string_view range1, std::string_view range2,
                                 std::string_view normalization, std::string_view count_type,
                                 std::string_view query_type) {
  if (range1.empty()) {
    return file_fetch_sum_all(f, normalization, count_type);
  }
  if (normalization != "NONE") {
    count_type = "float";
  }

  auto gi1 = query_type == "UCSC"
                 ? hictk::GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range1})
                 : hictk::GenomicInterval::parse_bed(f.chromosomes(), range1);
  auto gi2 = query_type == "UCSC"
                 ? hictk::GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range2})
                 : hictk::GenomicInterval::parse_bed(f.chromosomes(), range2);

  if (gi1 > gi2) {
    std::swap(gi1, gi2);
  }

  return std::visit(
      [&](const auto &ff) {
        auto sel =
            ff.fetch(gi1.chrom().name(), gi1.start(), gi1.end(), gi2.chrom().name(), gi2.start(),

                     gi2.end(), hictk::balancing::Method(normalization));

        if (count_type == "int") {
          return nb::cast(std::accumulate(
              sel.template begin<std::int32_t>(), sel.template end<std::int32_t>(), std::int64_t(0),
              [](const auto accumulator, const hictk::ThinPixel<std::int32_t> &p) {
                return accumulator + p.count;
              }));
        }
        return nb::cast(
            std::accumulate(sel.template begin<double>(), sel.template end<double>(), 0.0,
                            [](const auto accumulator, const hictk::ThinPixel<double> &p) {
                              return accumulator + p.count;
                            }));
      },
      f);
}

template <typename File>
inline std::int64_t file_fetch_nnz_all(const File &f) {
  if constexpr (std::is_same_v<File, hictk::cooler::File>) {
    return static_cast<std::int64_t>(f.nnz());
  }

  return std::visit(
      [&](const auto &ff) {
        auto sel = ff.fetch();
        return std::distance(sel.template begin<std::int32_t>(), sel.template end<std::int32_t>());
      },
      f);
}

template <typename File>
inline std::int64_t file_fetch_nnz(const File &f, std::string_view range1, std::string_view range2,
                                   std::string_view query_type) {
  if (range1.empty()) {
    return file_fetch_nnz_all(f);
  }
  auto gi1 = query_type == "UCSC"
                 ? hictk::GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range1})
                 : hictk::GenomicInterval::parse_bed(f.chromosomes(), range1);
  auto gi2 = query_type == "UCSC"
                 ? hictk::GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range2})
                 : hictk::GenomicInterval::parse_bed(f.chromosomes(), range2);

  if (gi1 > gi2) {
    std::swap(gi1, gi2);
  }

  return std::visit(
      [&](const auto &ff) {
        auto sel = f.fetch(gi1.chrom().name(), gi1.start(), gi1.end(), gi2.chrom().name(),
                           gi2.start(), gi2.end());
        return std::distance(sel.template begin<std::int32_t>(), sel.template end<std::int32_t>());
      },
      f);
}

}  // namespace hictkpy
