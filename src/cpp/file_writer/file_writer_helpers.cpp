// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/file_writer_helpers.hpp"

#include <arrow/array/array_base.h>
#include <arrow/array/builder_dict.h>
#include <arrow/array/builder_primitive.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <arrow/type_traits.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
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

namespace nb = nanobind;

namespace hictkpy {

template <typename Container>
[[nodiscard]] static bool str_contains(const Container& buff, std::string_view query) {
  const auto match =
      std::find_if(buff.begin(), buff.end(), [&](const auto& x) { return x == query; });
  return match != buff.end();
}

[[nodiscard]] static bool starts_with(std::string_view str, std::string_view prefix) noexcept {
  return str.find(prefix) == 0;
}

// clang-format off
using PixelDFColumn =
  std::variant<
    std::vector<std::string>,
    std::vector<std::int64_t>,
    std::vector<double>
>;
// clang-format on

[[nodiscard]] static bool issubdtype(const nb::module_& np, const nb::type_object& dtype1,
                                     const char* dtype2) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  return nb::cast<bool>(np.attr("issubdtype")(dtype1, np.attr(dtype2)));
}

[[nodiscard]] static bool is_integral(const nb::type_object& dtype) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  try {
    const auto np = import_module_checked("numpy");
    return issubdtype(np, dtype, "integer");
  } catch (...) {
    return nb::isinstance(nb::cast(std::int64_t{}), dtype);
  }
}

[[nodiscard]] static bool is_floating(const nb::type_object& dtype) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  try {
    const auto np = import_module_checked("numpy");
    return issubdtype(np, dtype, "floating");
  } catch (...) {
    return nb::isinstance(nb::cast(double{}), dtype);
  }
}

template <typename ArrowType>
[[nodiscard]] static std::shared_ptr<arrow::Array> make_array_of_numbers(
    const nb::object& iterable) {
  using N = typename arrow::TypeTraits<ArrowType>::CType;
  constexpr std::string_view type{std::is_floating_point_v<N> ? "number" : "integer"};

  try {
    arrow::NumericBuilder<ArrowType> builder;
    HICTKPY_GIL_SCOPED_ACQUIRE
    for (const auto& x : iterable) {
      try {
        const auto status = builder.Append(nb::cast<N>(x));
        if (!status.ok()) {
          const auto repr = nb::repr(x);
          throw std::runtime_error(
              fmt::format(FMT_STRING("failed to append {} to an array of {}: {}"), repr.c_str(),
                          type, status.message()));
        }
      } catch (const std::exception& e) {
        const auto repr = nb::repr(x);
        throw std::invalid_argument{
            fmt::format(FMT_STRING("failed to cast {} to a {}: {}"), repr.c_str(), type, e.what())};
      }
    }

    auto res = builder.Finish();
    if (!res.status().ok()) {
      throw std::runtime_error(res.status().message());
    }
    return res.MoveValueUnsafe();
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("failed to convert iterable to an array of {}: {}"), type, e.what()));
  }
}

[[nodiscard]] static std::shared_ptr<arrow::Array> make_array_of_integers(
    const nb::object& iterable) {
  return make_array_of_numbers<arrow::Int64Type>(iterable);
}

[[nodiscard]] static std::shared_ptr<arrow::Array> make_array_of_doubles(
    const nb::object& iterable) {
  return make_array_of_numbers<arrow::DoubleType>(iterable);
}

[[nodiscard]] static std::shared_ptr<arrow::Array> make_array_of_numbers(const nb::object& col) {
  try {
    HICTKPY_GIL_SCOPED_ACQUIRE
    std::optional<nb::type_object> type{};
    if (nb::hasattr(col, "dtype")) {
      type = nb::cast<nb::type_object>(col.attr("dtype"));
    } else if (nb::hasattr(col, "__getitem__")) {
      type = nb::cast<nb::type_object>(col.attr("__getitem__")(0).type());
    }

    if (!type.has_value()) {
      throw std::runtime_error("unable to cast object to an array of numbers: unknown object type");
    }
    if (is_integral(*type)) {
      return make_array_of_integers(col);
    }
    if (is_floating(*type)) {
      return make_array_of_doubles(col);
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to convert {} to int or float"), format_py_type(*type)));
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to cast object to an array of numbers: {}"), e.what()));
  } catch (...) {  // NOLINT
    throw std::runtime_error("unable to cast object to an array of numbers");
  }
}

[[nodiscard]] static std::shared_ptr<arrow::Array> make_array_of_strings(
    const nb::object& iterable) {
  try {
    arrow::StringDictionary32Builder builder;
    HICTKPY_GIL_SCOPED_ACQUIRE
    for (const auto& x : iterable) {
      try {
        const auto status = builder.Append(nb::cast<std::string>(x));
        if (!status.ok()) {
          const auto repr = nb::repr(x);
          throw std::runtime_error(
              fmt::format(FMT_STRING("failed to append {} to an array of strings: {}"),
                          repr.c_str(), status.message()));
        }
      } catch (const std::exception& e) {
        const auto repr = nb::repr(x);
        throw std::invalid_argument(
            fmt::format(FMT_STRING("failed to cast {} to a string: {}"), repr.c_str(), e.what()));
      }
    }

    auto res = builder.Finish();
    if (!res.status().ok()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to convert iterable to an array of strings: {}"),
                      res.status().message()));
    }
    return res.MoveValueUnsafe();
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to cast object to an array of strings: {}"), e.what()));
  } catch (...) {  // NOLINT
  }
  throw std::runtime_error("unable to cast object to an array of strings");
}

template <std::size_t N>
static void validate_columns(const std::vector<std::shared_ptr<arrow::Array>>& arrays,
                             const std::array<std::string_view, N>& column_names) {
  const auto size = arrays.front()->length();
  for (const auto& a : arrays) {
    if (a->length() != size) {
      std::vector<std::string> sizes;
      sizes.reserve(arrays.size());
      for (std::size_t i = 0; i < arrays.size(); ++i) {
        sizes.emplace_back(fmt::format(FMT_STRING("{}={}"), column_names[i], arrays[i]->length()));
      }
      throw std::invalid_argument(fmt::format(
          FMT_STRING("columns don't have the same lengths: [{}]"), fmt::join(sizes, ", ")));
    }
  }
}

template <std::size_t N>
[[nodiscard]] static std::vector<std::shared_ptr<arrow::Array>> preprocess_columns(
    const nb::dict& py_columns, const std::array<std::string_view, N>& column_names) {
  try {
    HICTKPY_GIL_SCOPED_ACQUIRE
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    arrays.reserve(py_columns.size());
    for (const auto& col_name : column_names) {
      try {
        auto col = py_columns.get(nb::str(col_name.data(), col_name.size()), nb::none());
        assert(!col.is_none());

        if (starts_with(col_name, "chrom")) {
          arrays.emplace_back(make_array_of_strings(col));
          continue;
        }
        if (col_name == "count") {
          arrays.emplace_back(make_array_of_numbers(col));
          continue;
        }
        arrays.emplace_back(make_array_of_integers(col));

      } catch (const std::exception& e) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("failed to process \"{}\" values: {}"), col_name, e.what()));
      } catch (...) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("failed to process \"{}\" values"), col_name));
      }
    }

    validate_columns(arrays, column_names);

    return arrays;
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to read pixels from dictionary: {}"), e.what()));
  } catch (...) {
    throw std::runtime_error("failed to read pixels from dictionary: unknown reason");
  }
}

[[nodiscard]] static std::vector<std::shared_ptr<arrow::Array>> preprocess_columns(
    PyArrowTable::Type table_type, const nb::dict& columns) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  if (table_type == PyArrowTable::Type::BG2) {
    return preprocess_columns(columns, bg2_columns);
  }

  if (table_type == PyArrowTable::Type::COO) {
    return preprocess_columns(columns, coo_columns);
  }
  unreachable_code();
}

[[nodiscard]] static PyArrowTable::Type infer_table_type(const nb::dict& columns) {
  std::uint8_t coo_cols_found{};
  std::uint8_t bed3_cols_found{};
  std::uint8_t bg2_cols_found{};

  HICTKPY_GIL_SCOPED_ACQUIRE
  for (const auto& py_col : columns.keys()) {
    const auto py_col_str = nb::str(py_col);
    const std::string_view col{py_col_str.c_str()};
    if (col == "count") {
      ++bg2_cols_found;
      ++coo_cols_found;
      continue;
    }
    if (str_contains(bg2_columns, col)) {
      ++bg2_cols_found;
      continue;
    }
    if (str_contains(coo_columns, col)) {
      ++coo_cols_found;
      continue;
    }
    if (str_contains(bed3_columns, col)) {
      ++bed3_cols_found;
    }
  }

  using T = PyArrowTable::Type;
  if (bg2_cols_found == bg2_columns.size()) {
    return T::BG2;
  }
  if (coo_cols_found == coo_columns.size()) {
    return T::COO;
  }
  if (bed3_cols_found == bed3_columns.size()) {
    return T::BED3;
  }
  return PyArrowTable::Type::UNKNOWN;
}

[[nodiscard]] static std::shared_ptr<arrow::Schema> make_schema(
    PyArrowTable::Type table_type, const std::shared_ptr<arrow::Array>& count_column) {
  arrow::FieldVector fields{};
  fields.reserve(table_type == PyArrowTable::Type::BG2 ? 7 : 3);  // NOLINT(*-avoid-magic-numbers)

  auto add_field = [&](const std::string& name, std::shared_ptr<arrow::DataType> type) {
    fields.emplace_back(std::make_shared<arrow::Field>(name, std::move(type)));
  };

  if (table_type == PyArrowTable::Type::COO) {
    add_field("bin1_id", arrow::int64());
    add_field("bin2_id", arrow::int64());
  } else if (table_type == PyArrowTable::Type::BG2) {
    add_field("chrom1", arrow::dictionary(arrow::int32(), arrow::utf8()));
    add_field("start1", arrow::int64());
    add_field("end1", arrow::int64());
    add_field("chrom2", arrow::dictionary(arrow::int32(), arrow::utf8()));
    add_field("start2", arrow::int64());
    add_field("end2", arrow::int64());
  }

  add_field("count", count_column->type());
  return std::make_shared<arrow::Schema>(std::move(fields));
}

namespace internal {

NumericDtype infer_count_type(const std::shared_ptr<arrow::Table>& df) {
  const auto type = infer_column_dtype(df, "count");
  return std::visit(
      [&df](const auto& x) -> NumericDtype {
        using T = remove_cvref_t<decltype(x)>;
        if constexpr (std::is_arithmetic_v<T>) {
          return x;
        } else {
          const auto col = df->GetColumnByName("count");
          if (!col) {
            throw std::runtime_error(
                "unable to infer dtype for column \"count\": column does not exist!");
          }
          throw std::runtime_error(fmt::format(
              FMT_STRING("unable to infer dtype for column \"count\": unable to map type "
                         "\"{}\" to a known numeric type"),
              col->type()->ToString()));
        }
      },
      type);
}

void raise_invalid_dict_format() {
  // clang-format off
  throw std::invalid_argument(
      "Dictionary does not contain columns in COO or BG2 format.\n"
      "Please make sure that the dictionary has the following keys:\n"
      "- COO: [bin1_id, bin2_id, count]\n"
      "- BG2: [chrom1, start1, end1, chrom2, start2, end2, count]\n"
      "And that values are iterable (e.g., list or numpy.array) with values of appropriate types:\n"
      "[chrom1, chrom2] -> string\n"
      "[bin1_id, bin2_id, start1, end1, start2, end2] -> int\n"
      "[count] -> int or float");
  // clang-format on
}

void raise_invalid_table_format() {
  // clang-format off
  throw std::invalid_argument(
      "DataFrame is not in COO or BG2 format.\n"
      "Please make sure that the DataFrame contains the following columns:\n"
      "- COO: [bin1_id, bin2_id, count]\n"
      "- BG2: [chrom1, start1, end1, chrom2, start2, end2, count]\n"
      "And that columns have appropriate dtypes:\n"
      "[chrom1, chrom2] -> string/categorical[string]\n"
      "[bin1_id, bin2_id, start1, end1, start2, end2] -> integral\n"
      "[count] -> numeric (excluding complex numbers)");
  // clang-format on
}

PyArrowTable make_table(const nb::dict& columns) {
  const auto table_type = infer_table_type(columns);
  if (table_type != PyArrowTable::Type::COO && table_type != PyArrowTable::Type::BG2) {
    raise_invalid_dict_format();
  }

  auto arrays = preprocess_columns(table_type, columns);
  auto schema = make_schema(table_type, arrays.back());

  return {arrow::Table::Make(std::move(schema), arrays), table_type};
}

}  // namespace internal

}  // namespace hictkpy
