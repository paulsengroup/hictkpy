// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/type.hpp"

#include <fmt/format.h>

#include <cstdint>
#include <exception>
#include <hictk/numeric_variant.hpp>
#include <string>
#include <string_view>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"

namespace hictkpy {

namespace nb = nanobind;

[[noreturn]] static void throw_exception(std::string_view dtype) {
  static const auto msg =
      fmt::format(FMT_STRING("Unable to map \"{}\" to a numeric C++ type.\n"
                             "Valid types are: uint, int, float, double, uint8, uint16, uint32, "
                             "uint64, int8, int16, int32, int64, float32, and float64."),
                  dtype);
  throw nb::type_error(msg.c_str());
}

[[nodiscard]] static nb::type_object get_dtype(const nb::module_& m, const char* name) {
  return m.attr(name);
}

[[nodiscard]] static hictk::internal::NumericVariant map_py_type_to_cpp_type(
    const nb::module_& m, const nb::type_object& dtype) {
  if (dtype.is(get_dtype(m, "uint8"))) {
    return std::uint8_t{};
  }
  if (dtype.is(get_dtype(m, "uint16"))) {
    return std::uint16_t{};
  }
  if (dtype.is(get_dtype(m, "uint32"))) {
    return std::uint32_t{};
  }
  if (dtype.is(get_dtype(m, "uint64"))) {
    return std::uint64_t{};
  }

  if (dtype.is(get_dtype(m, "int8"))) {
    return std::int8_t{};
  }
  if (dtype.is(get_dtype(m, "int16"))) {
    return std::int16_t{};
  }
  if (HICTKPY_LIKELY(dtype.is(get_dtype(m, "int32")))) {
    return std::int32_t{};
  }
  if (HICTKPY_LIKELY(dtype.is(get_dtype(m, "int64")))) {
    return std::int64_t{};
  }

  if (dtype.is(get_dtype(m, "float16"))) {
    return float{};
  }
  if (dtype.is(get_dtype(m, "float32"))) {
    return float{};
  }
  if (HICTKPY_LIKELY(dtype.is(get_dtype(m, "float64")))) {
    return double{};
  }

  const auto dtype_str = nb::str(dtype);
  throw_exception(dtype_str.c_str());
}

[[nodiscard]] static std::string nb_type_object_to_str(const nb::type_object& dtype) {
#ifdef _MSC_VER
  const auto dtype_mod = nb::str(nb::handle(dtype.attr("__module__")));
  const auto dtype_name = nb::str(nb::handle(dtype.attr("__name__")));
#else
  const auto dtype_mod = nb::str(dtype.attr("__module__"));
  const auto dtype_name = nb::str(dtype.attr("__name__"));
#endif

  if (std::string_view{dtype_mod.c_str()}.empty()) {
    return dtype_name.c_str();
  }
  return fmt::format(FMT_STRING("{}.{}"), dtype_mod.c_str(), dtype_name.c_str());
}

hictk::internal::NumericVariant map_py_type_to_cpp_type(const nb::type_object& dtype) {
  if (HICTKPY_LIKELY(dtype.is(nb::type<nb::int_>()))) {
    return std::int32_t{};
  }
  if (HICTKPY_LIKELY(dtype.is(nb::type<nb::float_>()))) {
    return double{};
  }

  auto np = import_module_checked("numpy");
  std::string_view exc_msg{};
  try {
    return map_py_type_to_cpp_type(np, dtype);
  } catch (const std::exception& e) {
    exc_msg = e.what();
  }

  // Sometimes checking Python's built-in types fails.
  // When this happens, fallback on comparisons based on type names
  try {
    return map_py_type_to_cpp_type(nb_type_object_to_str(dtype));
  } catch (const std::exception&) {
    throw nb::type_error(exc_msg.data());  // NOLINT(*-suspicious-stringview-data-usage)
  }
}

hictk::internal::NumericVariant map_py_type_to_cpp_type(std::string_view dtype) {
  // NOLINTBEGIN(*-avoid-magic-numbers)
  static_assert(sizeof(unsigned) == 4);
  static_assert(sizeof(int) == 4);
  static_assert(sizeof(float) == 4);
  static_assert(sizeof(double) == 8);
  // NOLINTEND(*-avoid-magic-numbers)

  if (const auto pos = dtype.rfind('.'); pos != std::string_view::npos) {
    try {
      return map_py_type_to_cpp_type(dtype.substr(pos + 1));
    } catch (const std::exception&) {
      throw_exception(dtype);
    }
  }

  if (dtype == "uint8") {
    return std::uint8_t{};
  }
  if (dtype == "uint16") {
    return std::uint16_t{};
  }
  if (dtype == "uint32" || dtype == "uint") {
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
  if (HICTKPY_LIKELY(dtype == "int32" || dtype == "int")) {
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
  if (HICTKPY_LIKELY(dtype == "float64" || dtype == "float" || dtype == "double")) {
    return double{};
  }

  throw_exception(dtype);
}

}  // namespace hictkpy
