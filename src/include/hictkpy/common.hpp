// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// clang-format off
#include "hictkpy/suppress_warnings.hpp"
HICTKPY_DISABLE_WARNING_PUSH
HICTKPY_DISABLE_WARNING_CAST_ALIGN
HICTKPY_DISABLE_WARNING_CXX98_COMPAT
HICTKPY_DISABLE_WARNING_OLD_STYLE_CAST
HICTKPY_DISABLE_WARNING_PEDANTIC
HICTKPY_DISABLE_WARNING_SHADOW
HICTKPY_DISABLE_WARNING_SIGN_CONVERSION
HICTKPY_DISABLE_WARNING_USELESS_CAST
#include <arrow/python/api.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
HICTKPY_DISABLE_WARNING_POP
// clang-format on

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
#include "hictk/transformers/join_genomic_coords.hpp"
#include "hictk/transformers/to_dataframe.hpp"

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
  using BufferT = nb::ndarray<nb::numpy, nb::shape<-1>, T>;
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
  inline void grow() { resize(static_cast<std::int64_t>(_buff.size()) * 2); }
  inline void shrink_to_fit() { resize(_size); }
  [[nodiscard]] auto operator()() -> BufferT {
    shrink_to_fit();
    return _buff;
  }
};

template <typename File>
inline nb::typed<nb::dict, std::string, std::uint32_t> get_chromosomes_from_file(
    const File &f, bool include_all = false) {
  nb::typed<nb::dict, std::string, std::uint32_t> py_chroms{};  // NOLINT
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

}  // namespace hictkpy
