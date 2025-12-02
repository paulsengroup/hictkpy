// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/type.hpp"

#include <fmt/format.h>

#include <cstdint>
#include <exception>
#include <string>
#include <string_view>

#include "hictkpy/common.hpp"
#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/variant.hpp"

namespace hictkpy {

namespace nb = nanobind;

[[nodiscard]] static std::string dtype_object_to_str(const nb::type_object& dtype) {
  try {
    HICTKPY_GIL_SCOPED_ACQUIRE
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
  } catch (const std::exception&) {
    return "unknown";
  }
}

[[noreturn]] static void throw_exception(std::string_view dtype, const char* msg = nullptr) {
  const auto s =
      fmt::format(FMT_STRING("Unable to map \"{}\" to a numeric C++ type.\n"
                             "Valid types are: uint, int, float, double, uint8, uint16, uint32, "
                             "uint64, int8, int16, int32, int64, float32, and float64.{}{}"),
                  dtype, msg ? "\n" : "", msg ? msg : "");

  HICTKPY_GIL_SCOPED_ACQUIRE
  throw nb::type_error(s.c_str());
}

[[noreturn]] static void throw_exception(const nb::type_object& dtype, const char* msg = nullptr) {
  throw_exception(dtype_object_to_str(dtype), msg);
}

[[nodiscard]] static bool issubdtype(const nb::module_& np, const nb::type_object& dtype1,
                                     const char* dtype2) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  return nb::cast<bool>(np.attr("issubdtype")(dtype1, np.attr(dtype2)));
}

NumericDtype map_py_numeric_to_cpp_type(const nb::type_object& dtype) {
  const auto np = import_module_checked("numpy");
  if (!issubdtype(np, dtype, "number")) {
    throw_exception(dtype, "Not a subdtype of numpy.number.");
  }

  if (issubdtype(np, dtype, "uint8")) {
    return std::uint8_t{};
  }
  if (issubdtype(np, dtype, "uint16")) {
    return std::uint16_t{};
  }
  if (issubdtype(np, dtype, "uint32")) {
    return std::uint32_t{};
  }
  if (issubdtype(np, dtype, "uint64")) {
    return std::uint64_t{};
  }

  if (issubdtype(np, dtype, "int8")) {
    return std::int8_t{};
  }
  if (issubdtype(np, dtype, "int16")) {
    return std::int16_t{};
  }
  if (HICTKPY_LIKELY(issubdtype(np, dtype, "int32"))) {
    return std::int32_t{};
  }
  if (HICTKPY_LIKELY(issubdtype(np, dtype, "int64"))) {
    return std::int64_t{};
  }

  if (issubdtype(np, dtype, "float32")) {
    return float{};
  }
  if (HICTKPY_LIKELY(issubdtype(np, dtype, "float64"))) {
    return double{};
  }

  throw_exception(dtype);
}

NumericDtype map_py_numeric_to_cpp_type(std::string_view dtype) {
  // NOLINTBEGIN(*-avoid-magic-numbers)
  static_assert(sizeof(unsigned) == 4);
  static_assert(sizeof(int) == 4);
  static_assert(sizeof(float) == 4);
  static_assert(sizeof(double) == 8);
  // NOLINTEND(*-avoid-magic-numbers)

  if (const auto pos = dtype.rfind('.'); pos != std::string_view::npos) {
    try {
      return map_py_numeric_to_cpp_type(dtype.substr(pos + 1));
    } catch (const std::exception& e) {
      throw_exception(dtype, e.what());
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

std::string format_py_type(const nb::handle& h) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  if (nb::hasattr(h, "__name__")) {
    return nb::cast<std::string>(h.attr("__name__"));
  }

  const auto type = nb::cast<nb::type_object>(h);
  return nb::cast<std::string>(type);
}

}  // namespace hictkpy
