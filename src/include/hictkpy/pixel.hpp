// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <hictk/bin.hpp>
#include <hictk/numeric_variant.hpp>
#include <hictk/pixel.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/type.hpp"

namespace hictkpy {

class Pixel {
  std::optional<hictk::PixelCoordinates> _coords{};
  std::int64_t _bin1_id{};
  std::int64_t _bin2_id{};
  std::variant<std::int64_t, double> _count{std::int64_t{}};

 public:
  Pixel() = default;
  template <typename N>
  Pixel(hictk::Pixel<N> p) noexcept  // NOLINT(*-explicit-conversions)
      : _coords(std::move(p.coords)),
        _bin1_id(static_cast<std::int64_t>(_coords->bin1.id())),
        _bin2_id(static_cast<std::int64_t>(_coords->bin2.id())),
        _count(cast_count(p.count)) {}

  template <typename N>
  Pixel(const hictk::ThinPixel<N> &p) noexcept  // NOLINT(*-explicit-conversions)
      : Pixel(static_cast<std::int64_t>(p.bin1_id), static_cast<std::int64_t>(p.bin2_id), p.count) {
  }

  template <typename N>
  Pixel(hictk::Bin bin1, hictk::Bin bin2, N count) noexcept
      : Pixel(hictk::Pixel{{std::move(bin1), std::move(bin2)}, count}) {}

  template <typename N>
  Pixel(std::int64_t bin1_id, std::int64_t bin2_id, N count_) noexcept
      : _bin1_id(bin1_id), _bin2_id(bin2_id), _count(cast_count(count_)) {}

  template <typename N>
  Pixel &operator=(hictk::Pixel<N> p) noexcept {
    _coords = std::move(p.coords);
    _bin1_id = static_cast<std::int64_t>(_coords->bin1.id());
    _bin2_id = static_cast<std::int64_t>(_coords->bin2.id());
    _count = cast_count(p.count);

    return *this;
  }

  template <typename N>
  Pixel &operator=(const hictk::ThinPixel<N> &p) noexcept {
    _coords = hictk::PixelCoordinates{};
    _bin1_id = static_cast<std::int64_t>(p.bin1_id);
    _bin2_id = static_cast<std::int64_t>(p.bin2_id);
    _count = cast_count(p.count);

    return *this;
  }

  [[nodiscard]] constexpr std::int64_t bin1_id() const noexcept { return _bin1_id; }
  [[nodiscard]] constexpr std::int64_t bin2_id() const noexcept { return _bin2_id; }
  [[nodiscard]] nanobind::object count() const {
    return std::visit([&](const auto n) { return nanobind::cast(n); }, _count);
  }

  [[nodiscard]] const hictk::PixelCoordinates &coords() const {
    if (HICTKPY_LIKELY(_coords.has_value())) {
      return *_coords;
    }

    throw nanobind::attribute_error(
        "Pixel does not have Bin with genomic coordinates associated with it. "
        "If you intend to access the genomic coordinates of Pixels, please make sure to call "
        "PixelSelector.fetch() with join=True.");
  }

  [[nodiscard]] const hictk::Bin &bin1() const { return coords().bin1; }
  [[nodiscard]] const hictk::Bin &bin2() const { return coords().bin2; }

  [[nodiscard]] std::string_view chrom1() const { return bin1().chrom().name(); }
  [[nodiscard]] std::int64_t start1() const { return static_cast<std::int64_t>(bin1().start()); }
  [[nodiscard]] std::int64_t end1() const { return static_cast<std::int64_t>(bin1().end()); }
  [[nodiscard]] std::string_view chrom2() const { return bin2().chrom().name(); }
  [[nodiscard]] std::int64_t start2() const { return static_cast<std::int64_t>(bin2().start()); }
  [[nodiscard]] std::int64_t end2() const { return static_cast<std::int64_t>(bin2().end()); }

  [[nodiscard]] std::string repr() const {
    return std::visit(
        [&](const auto &n) {
          if (_coords.has_value()) {
            return fmt::format(
                FMT_COMPILE("chrom1={}; start1={}; end1={}; chrom2={}; start2={}; end2={};"),
                _coords->bin1.chrom().name(), _coords->bin1.start(), _coords->bin1.end(),
                _coords->bin2.chrom().name(), _coords->bin2.start(), _coords->bin2.end(), n);
          }

          return fmt::format(FMT_COMPILE("bin1_id={}; bin2_id={}; count={};"), _bin1_id, _bin2_id,
                             n);
        },
        _count);
  }

  [[nodiscard]] std::string str() const {
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
    using Var = hictk::internal::NumericVariant;
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

  static void bind(nanobind::module_ &m) {
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

 private:
  template <typename N>
  [[nodiscard]] static constexpr std::variant<std::int64_t, double> cast_count(N n) noexcept {
    static_assert(std::is_arithmetic_v<N>);
    if constexpr (std::is_floating_point_v<N>) {
      return {conditional_static_cast<double>(n)};
    } else {
      return {conditional_static_cast<std::int64_t>(n)};
    }
  }
};

namespace internal {

template <typename T>
class ThreadSafeTypedPyList {
  static_assert(!std::is_same_v<nanobind::object, T>);
  std::optional<nanobind::list> _lst;

 public:
  ThreadSafeTypedPyList() : _lst(nanobind::list()) {}
  // NOLINTNEXTLINE(*-explicit-conversions)
  ThreadSafeTypedPyList(nanobind::list lst) noexcept : _lst(std::move(lst)) {}
  ThreadSafeTypedPyList(const ThreadSafeTypedPyList &) = default;
  ThreadSafeTypedPyList(ThreadSafeTypedPyList &&) noexcept = default;
  ~ThreadSafeTypedPyList() noexcept {
    HICTKPY_GIL_SCOPED_ACQUIRE
    _lst.reset();
  }

  ThreadSafeTypedPyList &operator=(const ThreadSafeTypedPyList &) = default;
  ThreadSafeTypedPyList &operator=(ThreadSafeTypedPyList &&) noexcept = default;

  [[nodiscard]] T at(std::size_t i) const {
    HICTKPY_GIL_SCOPED_ACQUIRE
    return nanobind::cast<T>((*_lst)[i]);  // NOLINT(*-unchecked-optional-access)
  }
};

template <typename N,
          typename NumpyArray = nanobind::ndarray<nanobind::numpy, nanobind::shape<-1>, N>>
[[nodiscard]] inline NumpyArray get_column_as_numpy(const nanobind::object &df,
                                                    std::string_view column) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  return nanobind::cast<NumpyArray>(df.attr("__getitem__")(column).attr("to_numpy")());
}

template <typename T>
[[nodiscard]] inline ThreadSafeTypedPyList<T> get_column_as_list(const nanobind::object &df,
                                                                 std::string_view column) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  return nanobind::cast<nanobind::list>(df.attr("__getitem__")(column).attr("to_list")());
}

inline void py_string_to_cpp(const nanobind::object &obj, std::string &buff) {
  buff.clear();
  HICTKPY_GIL_SCOPED_ACQUIRE
  auto sv = nanobind::cast<std::string_view>(obj);
  buff.assign(sv.data(), sv.size());
}

}  // namespace internal

template <typename N>
inline std::vector<hictk::ThinPixel<N>> coo_df_to_thin_pixels(const nanobind::object &df,
                                                              bool sort) {
  auto bin1_ids_np = internal::get_column_as_numpy<std::uint64_t>(df, "bin1_id");
  auto bin2_ids_np = internal::get_column_as_numpy<std::uint64_t>(df, "bin2_id");
  auto counts_np = internal::get_column_as_numpy<N>(df, "count");

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
                                                              const nanobind::object &df,
                                                              bool sort) {
  // TODO cast to pandas series with the appropriate string type?
  auto chrom1 = internal::get_column_as_list<std::string_view>(df, "chrom1");
  auto start1_np = internal::get_column_as_numpy<std::uint32_t>(df, "start1");
  auto end1_np = internal::get_column_as_numpy<std::uint32_t>(df, "end1");
  auto chrom2 = internal::get_column_as_list<std::string_view>(df, "chrom2");
  auto start2_np = internal::get_column_as_numpy<std::uint32_t>(df, "start2");
  auto end2_np = internal::get_column_as_numpy<std::uint32_t>(df, "end2");
  auto counts_np = internal::get_column_as_numpy<N>(df, "count");

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
          chrom1.at(i), start1(i), end1(i), chrom2.at(i), start2(i), end2(i), counts(i)));
    }

    auto bin1 = bin_table.at(reference.at(chrom1.at(i)), start1(i));
    auto bin2 = bin_table.at(reference.at(chrom2.at(i)), start2(i));

    if (bin_table.type() == hictk::BinTable::Type::fixed &&
        (end1(i) - start1(i) > bin_table.resolution() ||
         end2(i) - start2(i) > bin_table.resolution())) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Found an invalid pixel {} {} {} {} {} {} {}: pixel spans a "
                     "distance greater than the bin size"),
          chrom1.at(i), start1(i), end1(i), chrom2.at(i), start2(i), end2(i), counts(i)));
    }

    buffer[i] = hictk::ThinPixel<N>{bin1.id(), bin2.id(), counts(i)};
  }

  if (sort) {
    std::sort(buffer.begin(), buffer.end());
  }

  return buffer;
}

}  // namespace hictkpy
