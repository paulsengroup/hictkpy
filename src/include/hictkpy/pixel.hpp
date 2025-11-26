// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/array/array_base.h>
#include <arrow/array/array_binary.h>
#include <arrow/array/array_dict.h>
#include <arrow/array/array_primitive.h>
#include <arrow/compute/api_vector.h>
#include <arrow/compute/cast.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <arrow/type_fwd.h>
#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <hictk/bin.hpp>
#include <hictk/numeric_variant.hpp>
#include <hictk/pixel.hpp>
#include <limits>
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
#include "hictkpy/table.hpp"
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

template <typename... ChunkedArrays>
[[nodiscard]] inline auto normalize_non_uniform_column_types(
    const std::shared_ptr<arrow::DataType> &result_type, const ChunkedArrays &...arrays) {
  const auto first_array = [](auto first, auto &&...) { return first; }(arrays...);

  auto dtype_is_different = [&first_array](const auto &a) {
    return a->type()->id() != first_array->type()->id();
  };

  const bool casting_needed = (dtype_is_different(arrays) || ...);

  if (!casting_needed) {
    return std::make_tuple(arrays...);
  }

  auto cast_array = [&](const auto &array) {
    if (!dtype_is_different(array)) {
      return array;
    }
    SPDLOG_DEBUG(FMT_STRING("casting array from {} to {}..."), array->type()->ToString(),
                 result_type->ToString());
    auto res = arrow::compute::Cast(array, result_type);
    if (!res.ok()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to cast array of type {} to type {}: {}"),
                      array->type()->ToString(), result_type->ToString(), res.status().message()));
    }
    return res.ValueUnsafe().chunked_array();
  };

  return std::make_tuple(cast_array(arrays)...);
}

[[nodiscard]] inline std::shared_ptr<arrow::Table> ensure_table_has_uniform_chunks(
    std::shared_ptr<arrow::Table> table) {
  const auto &columns = table->columns();
  if (columns.size() < 2) {
    return table;
  }

  auto combine_table_chunks = [&]() {
    SPDLOG_DEBUG("found uneven chunks while converting arrow::Table to hictk::ThinPixels");
    auto res = table->CombineChunks();
    if (!res.ok()) {
      throw std::runtime_error(fmt::format(FMT_STRING("failed to combine arrow::Table chunks: {}"),
                                           res.status().message()));
    }
    table.reset();
    return res.MoveValueUnsafe();
  };

  auto check_uneven_chunks = [length = columns.front()->length()](const auto &chunk) {
    return chunk->length() != length;
  };

  auto num_chunks = columns.front()->num_chunks();
  auto table_has_uneven_chunks =
      std::any_of(columns.begin(), columns.end(),
                  [num_chunks](const auto &col) { return col->num_chunks() != num_chunks; });

  if (table_has_uneven_chunks) {
    return combine_table_chunks();
  }

  std::vector<std::shared_ptr<arrow::Array>> chunks(columns.size());
  for (std::int32_t i = 0; i < num_chunks; ++i) {
    chunks.clear();
    for (const auto &col : columns) {
      chunks.emplace_back(col->chunk(i));
    }

    if (std::any_of(chunks.begin(), chunks.end(), check_uneven_chunks)) {
      return combine_table_chunks();
    }
  }

  return table;
}

[[nodiscard]] inline auto numeric_array_static_pointer_cast_helper(
    const std::shared_ptr<arrow::DataType> &dtype) {
  // clang-format off
  using ArrowArray =
      std::variant<arrow::UInt8Array, arrow::UInt16Array, arrow::UInt32Array, arrow::UInt64Array,
                   arrow::Int8Array,  arrow::Int16Array,  arrow::Int32Array,  arrow::Int64Array,
                   arrow::FloatArray, arrow::DoubleArray>;
  // clang-format on

  switch (dtype->id()) {
    using T = arrow::Type::type;
    case T::UINT8:
      return ArrowArray{arrow::UInt8Array(nullptr)};
    case T::UINT16:
      return ArrowArray{arrow::UInt16Array(nullptr)};
    case T::UINT32:
      return ArrowArray{arrow::UInt32Array(nullptr)};
    case T::UINT64:
      return ArrowArray{arrow::UInt64Array(nullptr)};
    case T::INT8:
      return ArrowArray{arrow::Int8Array(nullptr)};
    case T::INT16:
      return ArrowArray{arrow::Int16Array(nullptr)};
    case T::INT32:
      return ArrowArray{arrow::Int32Array(nullptr)};
    case T::INT64:
      return ArrowArray{arrow::Int64Array(nullptr)};
    case T::FLOAT:
      return ArrowArray{arrow::FloatArray(nullptr)};
    case T::DOUBLE:
      return ArrowArray{arrow::DoubleArray(nullptr)};
    default:
      throw std::invalid_argument(
          fmt::format(FMT_STRING("{} is not a valid numeric dtype"), dtype->ToString()));
  }
}

[[nodiscard]] inline auto integer_array_static_pointer_cast_helper(
    const std::shared_ptr<arrow::DataType> &dtype) {
  // clang-format off
  using ArrowArray =
      std::variant<arrow::UInt8Array, arrow::UInt16Array, arrow::UInt32Array, arrow::UInt64Array,
                   arrow::Int8Array,  arrow::Int16Array,  arrow::Int32Array,  arrow::Int64Array>;
  // clang-format on

  return std::visit(
      [&]([[maybe_unused]] const auto &a) -> ArrowArray {
        using Array = remove_cvref_t<decltype(a)>;
        if constexpr (std::is_same_v<arrow::FloatArray, Array> ||
                      std::is_same_v<arrow::DoubleArray, Array>) {
          throw std::invalid_argument(
              fmt::format(FMT_STRING("{} is not a valid integral dtype"), dtype->ToString()));
        } else {
          return {Array(nullptr)};
        }
      },
      numeric_array_static_pointer_cast_helper(dtype));
}

[[nodiscard]] inline auto string_array_static_pointer_cast_helper(
    const std::shared_ptr<arrow::DataType> &dtype) {
  using ArrowArray =
      std::variant<arrow::StringArray, arrow::StringViewArray, arrow::LargeStringArray>;
  switch (dtype->id()) {
    using T = arrow::Type::type;
    case T::STRING:
      return ArrowArray{arrow::StringArray(nullptr)};
    case T::STRING_VIEW:
      return ArrowArray{arrow::StringViewArray(nullptr)};
    case T::LARGE_STRING:
      return ArrowArray{arrow::LargeStringArray(nullptr)};
    default:
      throw std::invalid_argument(
          fmt::format(FMT_STRING("{} is not a valid string dtype"), dtype->ToString()));
  }
}

namespace coo {
template <typename ArrowBin, typename ArrowCount, typename Count>
static void process_chunk(const std::shared_ptr<arrow::NumericArray<ArrowBin>> &bin1_ids_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &bin2_ids_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  const auto size = bin1_ids_chunk->data()->length;
  assert(bin2_ids_chunk->data()->length == size);
  assert(counts_chunk->data()->length == size);

  if (size == 0) {
    return;
  }

  for (std::int64_t i = 0; i < size; ++i) {
    const auto b1 = bin1_ids_chunk->Value(i);
    const auto b2 = bin2_ids_chunk->Value(i);
    const auto count = counts_chunk->Value(i);

    if constexpr (std::is_signed_v<decltype(b1)>) {
      if (b1 < 0) {
        throw std::runtime_error("found negative value in bin1_id column");
      }
      if (b2 < 0) {
        throw std::runtime_error("found negative value in bin2_id column");
      }
    }
    buff.emplace_back(hictk::ThinPixel<Count>{conditional_static_cast<std::uint64_t>(b1),
                                              conditional_static_cast<std::uint64_t>(b2),
                                              conditional_static_cast<Count>(count)});
  }
}

template <typename ArrowCount, typename Count>
static void process_chunk(const std::shared_ptr<arrow::Array> &bin1_ids_chunk,
                          const std::shared_ptr<arrow::Array> &bin2_ids_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  assert(bin1_ids_chunk->type()->id() == bin2_ids_chunk->type()->id());

  const auto array_type_var = [&]() {
    try {
      return integer_array_static_pointer_cast_helper(bin1_ids_chunk->type());
    } catch (const std::exception &e) {
      throw std::invalid_argument(
          fmt::format(FMT_STRING("failed to infer dtype for bin{{1,2}}_id columns: {}"), e.what()));
    } catch (...) {
      throw std::invalid_argument("failed to infer dtype for bin{{1,2}}_id columns: unknown error");
    }
  }();

  return std::visit(
      [&]([[maybe_unused]] const auto &array) {
        using T = remove_cvref_t<decltype(array)>;
        process_chunk(std::static_pointer_cast<T>(bin1_ids_chunk),
                      std::static_pointer_cast<T>(bin2_ids_chunk), counts_chunk, buff);
      },
      array_type_var);
}

template <typename Count>
static void process_chunk(const std::shared_ptr<arrow::Array> &bin1_ids_chunk,
                          const std::shared_ptr<arrow::Array> &bin2_ids_chunk,
                          const std::shared_ptr<arrow::Array> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  const auto array_type_var = [&]() {
    try {
      return numeric_array_static_pointer_cast_helper(counts_chunk->type());
    } catch (const std::exception &e) {
      throw std::invalid_argument(
          fmt::format(FMT_STRING("failed to infer dtype for count column: {}"), e.what()));
    } catch (...) {
      throw std::invalid_argument("failed to infer dtype for count column: unknown error");
    }
  }();

  return std::visit(
      [&]([[maybe_unused]] const auto &array) {
        using T = remove_cvref_t<decltype(array)>;
        process_chunk(bin1_ids_chunk, bin2_ids_chunk, std::static_pointer_cast<T>(counts_chunk),
                      buff);
      },
      array_type_var);
}

template <typename N>
inline std::vector<hictk::ThinPixel<N>> convert_table_thin_pixels(std::shared_ptr<arrow::Table> df,
                                                                  bool sort) {
  try {
    df = ensure_table_has_uniform_chunks(df);

    // we assume that the array types have already been validated
    auto [bin1_ids, bin2_ids] = normalize_non_uniform_column_types(
        arrow::int64(), df->GetColumnByName("bin1_id"), df->GetColumnByName("bin2_id"));
    auto counts = df->GetColumnByName("count");

    std::vector<hictk::ThinPixel<N>> buffer;
    buffer.reserve(static_cast<std::size_t>(df->num_rows()));

    const auto num_chunks = counts->num_chunks();
    for (std::int32_t i = 0; i < num_chunks; ++i) {
      const auto chunk1 = bin1_ids->chunk(i);
      const auto chunk2 = bin2_ids->chunk(i);
      const auto chunk3 = counts->chunk(i);
      process_chunk(chunk1, chunk2, chunk3, buffer);
    }

    assert(buffer.size() == static_cast<std::size_t>(df->num_rows()));
    df.reset();

    if (sort) {
      std::sort(buffer.begin(), buffer.end());
    }

    return buffer;
  } catch (const std::invalid_argument &e) {
    throw std::invalid_argument(
        fmt::format(FMT_STRING("failed to convert DataFrame to a COO pixels: {}"), e.what()));

  } catch (const std::runtime_error &e) {
    throw std::invalid_argument(
        fmt::format(FMT_STRING("failed to convert DataFrame to a COO pixels: {}"), e.what()));
  } catch (...) {
    throw std::invalid_argument("failed to convert DataFrame to a COO pixels: unknown error");
  }
}

}  // namespace coo

namespace bg2 {

template <typename I, typename ArrowStringLikeArray>
class ChromosomeArray {
  std::shared_ptr<arrow::NumericArray<I>> _indices{};
  std::shared_ptr<ArrowStringLikeArray> _names{};
  const hictk::Reference *_chroms{};

 public:
  ChromosomeArray() = default;
  ChromosomeArray(const hictk::BinTable &bins, std::shared_ptr<arrow::NumericArray<I>> indices,
                  std::shared_ptr<ArrowStringLikeArray> names)
      : _indices(std::move(indices)), _names(std::move(names)), _chroms(&bins.chromosomes()) {}

  [[nodiscard]] std::size_t size() const noexcept {
    if (!_indices) {
      return 0;
    }
    return static_cast<std::size_t>(_indices->length());
  }

  [[nodiscard]] std::string_view get_chrom_name(std::int64_t i) const {
    if (!_indices) {
      throw std::out_of_range("out of range: ChromosomeArray is empty");
    }

    assert(i < size());

    const auto idx = conditional_static_cast<std::int64_t>(_indices->Value(i));
    return _names->GetView(idx);
  }

  [[nodiscard]] std::uint32_t get_chrom_id(std::int64_t i) const {
    if (!_chroms) {
      return std::numeric_limits<std::uint32_t>::max();
    }
    return _chroms->at(get_chrom_name(i)).id();
  }
};

template <typename ArrowBin, typename ArrowCount, typename Count, typename ArrowChromIdx1,
          typename ArrowChromStr1, typename ArrowChromIdx2, typename ArrowChromStr2>
static void process_chunk(const hictk::BinTable &bins,
                          const ChromosomeArray<ArrowChromIdx1, ArrowChromStr1> &chrom1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end1_chunk,
                          const ChromosomeArray<ArrowChromIdx2, ArrowChromStr2> &chrom2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  const auto size = static_cast<std::int64_t>(chrom1_chunk.size());
  assert(start1_chunk->length() == size);
  assert(end1_chunk->length() == size);
  assert(chrom2_chunk.size() == size);
  assert(start2_chunk->length() == size);
  assert(end2_chunk->length() == size);
  assert(counts_chunk->length() == size);

  auto get_bin_checked = [&](std::uint_fast8_t idx, auto chrom_id, auto start, auto end) {
    auto bin = bins.at(chrom_id, start);
    if (bin.end() == end) {
      return bin;
    }
    const auto res = bins.resolution();
    throw std::runtime_error(fmt::format(
        FMT_STRING("invalid end{}, expected {}, found {}, (start{}={}; bin_size={})"), idx,
        bin.end(), end, idx, bin.start(), res == 0 ? "variable" : fmt::to_string(res)));
  };

  try {
    for (std::int64_t i = 0; i < size; ++i) {
      const auto chrom1_id = chrom1_chunk.get_chrom_id(i);
      const auto start1 = start1_chunk->Value(i);
      const auto end1 = end1_chunk->Value(i);
      const auto chrom2_id = chrom2_chunk.get_chrom_id(i);
      const auto start2 = start2_chunk->Value(i);
      const auto end2 = end2_chunk->Value(i);
      const auto count = counts_chunk->Value(i);
      try {
        if (start1 < 0 || end1 < 0 || start2 < 0 || end2 < 0) {
          throw std::runtime_error("genomic coordinates cannot be negative");
        }

        if (end1 < start1 || end2 < start2) {
          throw std::runtime_error(
              "end position of a bin cannot be smaller than its start position");
        }

        const hictk::Pixel p{
            get_bin_checked(1, chrom1_id, conditional_static_cast<std::uint32_t>(start1),
                            conditional_static_cast<std::uint32_t>(end1)),
            get_bin_checked(2, chrom2_id, conditional_static_cast<std::uint32_t>(start2),
                            conditional_static_cast<std::uint32_t>(end2)),
            conditional_static_cast<Count>(count)};
        buff.emplace_back(p.to_thin());
      } catch (const std::exception &e) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("failed to map {}:{}-{}; {}:{}-{} to a valid pixel: {}"),
                        chrom1_chunk.get_chrom_name(i), start1, end1,
                        chrom2_chunk.get_chrom_name(i), start1, end1, e.what()));
      } catch (...) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("failed to map {}:{}-{}; {}:{}-{} to a valid pixel: unknown error"),
            chrom1_chunk.get_chrom_name(i), start1, end1, chrom2_chunk.get_chrom_name(i), start1,
            end1));
      }
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to process BG2 pixel: {}"), e.what()));
  } catch (...) {
    throw std::runtime_error("failed to process BG2 pixel: unknown error");
  }
}

template <typename ArrowBin, typename ArrowCount, typename Count, typename ArrowChromIdx1,
          typename ArrowChromStr1>
static void process_chunk(const hictk::BinTable &bins,
                          const ChromosomeArray<ArrowChromIdx1, ArrowChromStr1> &chrom1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end1_chunk,
                          const std::shared_ptr<arrow::DictionaryArray> &chrom2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  std::visit(
      [&]([[maybe_unused]] const auto &array1) {
        using T1 = remove_cvref_t<decltype(array1)>;
        std::visit(
            [&]([[maybe_unused]] const auto &array2) {
              using T2 = remove_cvref_t<decltype(array2)>;
              process_chunk(
                  bins, chrom1_chunk, start1_chunk, end1_chunk,
                  ChromosomeArray(bins, std::static_pointer_cast<T1>(chrom2_chunk->indices()),
                                  std::static_pointer_cast<T2>(chrom2_chunk->dictionary())),
                  start2_chunk, end2_chunk, counts_chunk, buff);
            },
            string_array_static_pointer_cast_helper(chrom2_chunk->dictionary()->type()));
      },
      integer_array_static_pointer_cast_helper(chrom2_chunk->indices()->type()));
}

template <typename ArrowBin, typename ArrowCount, typename Count>
static void process_chunk(const hictk::BinTable &bins,
                          const std::shared_ptr<arrow::DictionaryArray> &chrom1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end1_chunk,
                          const std::shared_ptr<arrow::DictionaryArray> &chrom2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  std::visit(
      [&]([[maybe_unused]] const auto &array1) {
        using T1 = remove_cvref_t<decltype(array1)>;
        std::visit(
            [&]([[maybe_unused]] const auto &array2) {
              using T2 = remove_cvref_t<decltype(array2)>;
              process_chunk(
                  bins,
                  ChromosomeArray(bins, std::static_pointer_cast<T1>(chrom1_chunk->indices()),
                                  std::static_pointer_cast<T2>(chrom1_chunk->dictionary())),
                  start1_chunk, end1_chunk, chrom2_chunk, start2_chunk, end2_chunk, counts_chunk,
                  buff);
            },
            string_array_static_pointer_cast_helper(chrom1_chunk->dictionary()->type()));
      },
      integer_array_static_pointer_cast_helper(chrom1_chunk->indices()->type()));
}

template <typename ArrowBin, typename ArrowCount, typename Count>
static void process_chunk(const hictk::BinTable &bins,
                          const std::shared_ptr<arrow::Array> &chrom1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start1_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end1_chunk,
                          const std::shared_ptr<arrow::Array> &chrom2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &start2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowBin>> &end2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  auto encode = [](std::string_view name, const std::shared_ptr<arrow::Array> &a) {
    // this should be fine, as calling dictionary encoded on an already encoded array is a no-op
    auto res = arrow::compute::DictionaryEncode(a);
    if (!res.ok()) {
      throw std::invalid_argument(
          fmt::format(FMT_STRING("failed to dictionary encode values in {} column: {}"), name,
                      res.status().message()));
    }

    return res.MoveValueUnsafe().array_as<arrow::DictionaryArray>();
  };

  // clang-format off
  process_chunk(bins,
    encode("chrom1", chrom1_chunk), start1_chunk, end1_chunk,
    encode("chrom2", chrom2_chunk), start2_chunk, end2_chunk,
    counts_chunk,
    buff
  );
  // clang-format on
}

template <typename ArrowCount, typename Count>
static void process_chunk(const hictk::BinTable &bins,
                          const std::shared_ptr<arrow::Array> &chrom1_chunk,
                          const std::shared_ptr<arrow::Array> &start1_chunk,
                          const std::shared_ptr<arrow::Array> &end1_chunk,
                          const std::shared_ptr<arrow::Array> &chrom2_chunk,
                          const std::shared_ptr<arrow::Array> &start2_chunk,
                          const std::shared_ptr<arrow::Array> &end2_chunk,
                          const std::shared_ptr<arrow::NumericArray<ArrowCount>> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  assert(start1_chunk->type()->id() == end1_chunk->type()->id());
  assert(start1_chunk->type()->id() == start2_chunk->type()->id());
  assert(start1_chunk->type()->id() == end2_chunk->type()->id());

  const auto array_type_var = [&]() {
    try {
      return numeric_array_static_pointer_cast_helper(start1_chunk->type());
    } catch (const std::exception &e) {
      throw std::invalid_argument(fmt::format(
          FMT_STRING("failed to infer dtype for start{{1,2}} and end{{1,2}} columns: {}"),
          e.what()));
    } catch (...) {
      throw std::invalid_argument(
          "failed to infer dtype for start{{1,2}} and end{{1,2}} columns: unknown error");
    }
  }();

  return std::visit(
      [&]([[maybe_unused]] const auto &array) {
        using T = remove_cvref_t<decltype(array)>;
        // clang-format off
        process_chunk(
            bins,
            chrom1_chunk,
            std::static_pointer_cast<T>(start1_chunk),
            std::static_pointer_cast<T>(end1_chunk),
            chrom2_chunk,
            std::static_pointer_cast<T>(start2_chunk),
            std::static_pointer_cast<T>(end2_chunk),
            counts_chunk,
            buff
        );
        // clang-format on
      },
      array_type_var);
}

template <typename Count>
static void process_chunk(const hictk::BinTable &bins,
                          const std::shared_ptr<arrow::Array> &chrom1_chunk,
                          const std::shared_ptr<arrow::Array> &start1_chunk,
                          const std::shared_ptr<arrow::Array> &end1_chunk,
                          const std::shared_ptr<arrow::Array> &chrom2_chunk,
                          const std::shared_ptr<arrow::Array> &start2_chunk,
                          const std::shared_ptr<arrow::Array> &end2_chunk,
                          const std::shared_ptr<arrow::Array> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  const auto array_type_var = [&]() {
    try {
      return numeric_array_static_pointer_cast_helper(counts_chunk->type());
    } catch (const std::exception &e) {
      throw std::invalid_argument(
          fmt::format(FMT_STRING("failed to infer dtype for count column: {}"), e.what()));
    } catch (...) {
      throw std::invalid_argument("failed to infer dtype for count column: unknown error");
    }
  }();

  return std::visit(
      [&]([[maybe_unused]] const auto &array) {
        using T = remove_cvref_t<decltype(array)>;
        process_chunk(bins, chrom1_chunk, start1_chunk, end1_chunk, chrom2_chunk, start2_chunk,
                      end2_chunk, std::static_pointer_cast<T>(counts_chunk), buff);
      },
      array_type_var);
}

template <typename N>
inline std::vector<hictk::ThinPixel<N>> convert_table_thin_pixels(const hictk::BinTable &bins,
                                                                  std::shared_ptr<arrow::Table> df,
                                                                  bool sort) {
  try {
    df = ensure_table_has_uniform_chunks(df);
    auto chrom1 = df->GetColumnByName("chrom1");
    auto chrom2 = df->GetColumnByName("chrom2");
    auto counts = df->GetColumnByName("count");
    auto [start1, end1, start2, end2] = normalize_non_uniform_column_types(
        arrow::int64(), df->GetColumnByName("start1"), df->GetColumnByName("end1"),
        df->GetColumnByName("start2"), df->GetColumnByName("end2"));

    std::vector<hictk::ThinPixel<N>> buffer;
    buffer.reserve(static_cast<std::size_t>(df->num_rows()));

    const auto num_chunks = counts->num_chunks();
    for (std::int32_t i = 0; i < num_chunks; ++i) {
      auto chunk1 = chrom1->chunk(i);
      auto chunk2 = start1->chunk(i);
      auto chunk3 = end1->chunk(i);
      auto chunk4 = chrom2->chunk(i);
      auto chunk5 = start2->chunk(i);
      auto chunk6 = end2->chunk(i);
      auto chunk7 = counts->chunk(i);
      process_chunk(bins, chunk1, chunk2, chunk3, chunk4, chunk5, chunk6, chunk7, buffer);
    }

    assert(buffer.size() == static_cast<std::size_t>(df->num_rows()));
    df.reset();

    if (sort) {
      std::sort(buffer.begin(), buffer.end());
    }

    return buffer;
  } catch (const std::invalid_argument &e) {
    throw std::invalid_argument(
        fmt::format(FMT_STRING("failed to convert DataFrame to a COO pixels: {}"), e.what()));

  } catch (const std::runtime_error &e) {
    throw std::invalid_argument(
        fmt::format(FMT_STRING("failed to convert DataFrame to a COO pixels: {}"), e.what()));
  } catch (...) {
    throw std::invalid_argument("failed to convert DataFrame to a COO pixels: unknown error");
  }
}

}  // namespace bg2

template <typename N>
inline std::vector<hictk::ThinPixel<N>> convert_table_to_thin_pixels(
    const hictk::BinTable &bin_table, const PyArrowTable &df, bool sort) {
  switch (df.type()) {
    case PyArrowTable::Type::COO:
      return coo::convert_table_thin_pixels<N>(df, sort);
    case PyArrowTable::Type::BG2:
      return bg2::convert_table_thin_pixels<N>(bin_table, df, sort);
    default:
      unreachable_code();
  }
}

}  // namespace hictkpy
