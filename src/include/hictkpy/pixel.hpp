// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/array/array_base.h>
#include <arrow/array/array_primitive.h>
#include <arrow/compute/cast.h>
#include <arrow/table.h>
#include <arrow/type.h>
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

  using T = arrow::Type::type;
  switch (bin1_ids_chunk->type()->id()) {
    case T::UINT8: {
      process_chunk(std::static_pointer_cast<arrow::UInt8Array>(bin1_ids_chunk),
                    std::static_pointer_cast<arrow::UInt8Array>(bin2_ids_chunk), counts_chunk,
                    buff);
      return;
    }
    case T::UINT16: {
      process_chunk(std::static_pointer_cast<arrow::UInt16Array>(bin1_ids_chunk),
                    std::static_pointer_cast<arrow::UInt16Array>(bin2_ids_chunk), counts_chunk,
                    buff);
      return;
    }
    case T::UINT32: {
      process_chunk(std::static_pointer_cast<arrow::UInt32Array>(bin1_ids_chunk),
                    std::static_pointer_cast<arrow::UInt32Array>(bin2_ids_chunk), counts_chunk,
                    buff);
      return;
    }
    case T::UINT64: {
      process_chunk(std::static_pointer_cast<arrow::UInt64Array>(bin1_ids_chunk),
                    std::static_pointer_cast<arrow::UInt64Array>(bin2_ids_chunk), counts_chunk,
                    buff);
      return;
    }
    case T::INT8: {
      process_chunk(std::static_pointer_cast<arrow::Int8Array>(bin1_ids_chunk),
                    std::static_pointer_cast<arrow::Int8Array>(bin2_ids_chunk), counts_chunk, buff);
      return;
    }
    case T::INT16: {
      process_chunk(std::static_pointer_cast<arrow::Int16Array>(bin1_ids_chunk),
                    std::static_pointer_cast<arrow::Int16Array>(bin2_ids_chunk), counts_chunk,
                    buff);
      return;
    }
    case T::INT32: {
      process_chunk(std::static_pointer_cast<arrow::Int32Array>(bin1_ids_chunk),
                    std::static_pointer_cast<arrow::Int32Array>(bin2_ids_chunk), counts_chunk,
                    buff);
      return;
    }
    case T::INT64: {
      process_chunk(std::static_pointer_cast<arrow::Int64Array>(bin1_ids_chunk),
                    std::static_pointer_cast<arrow::Int64Array>(bin2_ids_chunk), counts_chunk,
                    buff);
      return;
    }
    default:
      throw std::invalid_argument("failed to infer dtype for bin{1,2}_id columns: unknown error");
  }
}

template <typename Count>
static void process_chunk(const std::shared_ptr<arrow::Array> &bin1_ids_chunk,
                          const std::shared_ptr<arrow::Array> &bin2_ids_chunk,
                          const std::shared_ptr<arrow::Array> &counts_chunk,
                          std::vector<hictk::ThinPixel<Count>> &buff) {
  using T = arrow::Type::type;
  switch (counts_chunk->type()->id()) {
    case T::UINT8: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::UInt8Array>(counts_chunk), buff);
      return;
    }
    case T::UINT16: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::UInt16Array>(counts_chunk), buff);
      return;
    }
    case T::UINT32: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::UInt32Array>(counts_chunk), buff);
      return;
    }
    case T::UINT64: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::UInt64Array>(counts_chunk), buff);
      return;
    }
    case T::INT8: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::Int8Array>(counts_chunk), buff);
      return;
    }
    case T::INT16: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::Int16Array>(counts_chunk), buff);
      return;
    }
    case T::INT32: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::Int32Array>(counts_chunk), buff);
      return;
    }
    case T::INT64: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::Int64Array>(counts_chunk), buff);
      return;
    }
    case T::FLOAT: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::FloatArray>(counts_chunk), buff);
      return;
    }
    case T::DOUBLE: {
      process_chunk(bin1_ids_chunk, bin2_ids_chunk,
                    std::static_pointer_cast<arrow::DoubleArray>(counts_chunk), buff);
      return;
    }
    default:
      throw std::invalid_argument("failed to infer dtype for count column: unknown error");
  }
}

template <typename N>
inline std::vector<hictk::ThinPixel<N>> convert_table_thin_pixels(std::shared_ptr<arrow::Table> df,
                                                                  bool sort) {
  try {
    // we assume that the array types have already been validated
    auto bin1_ids = df->GetColumnByName("bin1_id");
    auto bin2_ids = df->GetColumnByName("bin2_id");
    auto counts = df->GetColumnByName("count");

    if (*bin1_ids->type() != *bin2_ids->type()) {
      if (*bin1_ids->type() != *arrow::int64()) {
        SPDLOG_DEBUG(FMT_STRING("casting bin1_id from {} to {}..."), bin1_ids->type()->ToString(),
                     arrow::int64()->ToString());
        auto res = arrow::compute::Cast(bin1_ids, arrow::int64());
        if (!res.ok()) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("failed to cast array of type {} to type {}: {}"),
              bin1_ids->type()->ToString(), arrow::int64()->ToString(), res.status().message()));
        }
        bin1_ids = res.ValueUnsafe().chunked_array();
      }

      if (*bin2_ids->type() != *arrow::int64()) {
        SPDLOG_DEBUG(FMT_STRING("casting bin2_id from {} to {}..."), bin2_ids->type()->ToString(),
                     arrow::int64()->ToString());
        auto res = arrow::compute::Cast(bin2_ids, arrow::int64());
        if (!res.ok()) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("failed to cast array of type {} to type {}: {}"),
              bin2_ids->type()->ToString(), arrow::int64()->ToString(), res.status().message()));
        }
        bin2_ids = res.ValueUnsafe().chunked_array();
      }
    }

    // ensure columns have the same number of chunks
    bool uneven_chunks = false;
    if (bin1_ids->num_chunks() != bin2_ids->num_chunks() ||
        bin1_ids->num_chunks() != counts->num_chunks()) {
      uneven_chunks = true;
    }

    const auto num_chunks = bin1_ids->num_chunks();
    using I = remove_cvref_t<decltype(num_chunks)>;
    if (!uneven_chunks) {
      // ensure chunks have the same size
      for (I i = 0; i < num_chunks; ++i) {
        const auto chunk1 = bin1_ids->chunk(i);
        const auto chunk2 = bin2_ids->chunk(i);
        const auto chunk3 = counts->chunk(i);

        if (chunk1->length() != chunk2->length() || chunk1->length() != chunk3->length()) {
          uneven_chunks = true;
          break;
        }
      }
    }

    if (uneven_chunks) {
      SPDLOG_DEBUG("found uneven chunks while converting arrow::Table to hictk::ThinPixels");
      auto res = df->CombineChunks();
      if (!res.ok()) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("failed to combine arrow::Table chunks: {}"), res.status().message()));
      }
      df = res.MoveValueUnsafe();
    }

    std::vector<hictk::ThinPixel<N>> buffer;
    buffer.reserve(static_cast<std::size_t>(df->num_rows()));

    for (I i = 0; i < num_chunks; ++i) {
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
