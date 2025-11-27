// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Python.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstdio>
#include <exception>
#include <string>
#include <utility>
#include <vector>

#include "hictkpy/locking.hpp"

namespace hictkpy {

[[nodiscard]] inline nanobind::module_ import_module_checked(const std::string& module_name) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  try {
    return nanobind::module_::import_(module_name.c_str());
  } catch (nanobind::python_error& e) {
    // NOLINTNEXTLINE(*-vararg)
    nanobind::raise_from(e, PyExc_ModuleNotFoundError,
                         "To enable %s support, please install %s with: pip install 'hictkpy[%s]'\n"
                         "Alternatively, you can install hictkpy with all its dependencies by "
                         "running: pip install 'hictkpy[all]'",
                         module_name.c_str(), module_name.c_str(), module_name.c_str());
  }
}

// NOLINTNEXTLINE(*-avoid-magic-numbers)
[[nodiscard]] inline nanobind::module_ import_pyarrow_checked(int min_version_major = 16,
                                                              int min_version_minor = 0,
                                                              int min_version_patch = 0) {
  assert(min_version_major >= 0);
  assert(min_version_minor >= 0);
  assert(min_version_patch >= 0);

  [[maybe_unused]] const GilScopedAcquire gil{true};

  static bool version_ok{false};
  auto pa = import_module_checked("pyarrow");

  if (version_ok) {
    return pa;
  }

  static std::string error_msg{};
  error_msg.clear();
  try {
    auto metadata = nanobind::module_::import_("importlib.metadata");

    const auto version = nanobind::cast<std::vector<std::string>>(
        metadata.attr("version")("pyarrow").attr("split")("."));
    if (version.size() < 3) {
      throw nanobind::import_error(
          "unable to detect pyarrow version: assuming pyarrow's version is not compatible: please "
          "install a compatible version of pyarrow with: pip install 'hictkpy[pyarrow]'");
    }

    const auto major_version_found = std::stoi(version[0]);
    const auto minor_version_found = std::stoi(version[1]);
    const auto patch_version_found = std::stoi(version[2]);

    version_ok = major_version_found >= min_version_major;
    version_ok |=
        major_version_found == min_version_major && minor_version_found >= min_version_minor;
    version_ok |= major_version_found == min_version_major &&
                  minor_version_found == min_version_minor &&
                  patch_version_found >= min_version_patch;

    if (!version_ok) {
      error_msg = fmt::format(
          FMT_STRING("pyarrow {}.{}.{} is too old to be used with hictkpy: please "
                     "install a compatible version with: pip install 'hictkpy[pyarrow]'"),
          major_version_found, minor_version_found, patch_version_found);
    }
  } catch (const nanobind::builtin_exception&) {
    throw;
  } catch (const std::exception& e) {
    error_msg = fmt::format(
        FMT_STRING(
            "Unable to parse pyarrow version: {}.\n"
            "Assuming pyarrow's version is not compatible.\n"
            "Please install a compatible version of pyarrow with: pip install 'hictkpy[pyarrow]'"),
        e.what());
  } catch (...) {
    error_msg =
        "Unable to parse pyarrow version.\n"
        "Assuming pyarrow's version is not compatible.\n"
        "Please install a compatible version of pyarrow with: pip install 'hictkpy[pyarrow]'";
  }

  if (!error_msg.empty()) {
    throw nanobind::import_error(error_msg.c_str());
  }

  assert(version_ok);
  return pa;
}

inline void check_module_is_importable(const std::string& module_name) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  std::ignore = import_module_checked(module_name);
}

// NOLINTNEXTLINE(*-avoid-magic-numbers)
inline void check_pyarrow_is_importable(int min_version_major = 16, int min_version_minor = 0,
                                        int min_version_patch = 0) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  std::ignore = import_pyarrow_checked(min_version_major, min_version_minor, min_version_patch);
}

template <typename... T>
inline void println_stderr_noexcept(fmt::format_string<T...> fmt, T&&... args) noexcept {
  try {
    fmt::println(stderr, fmt, std::forward<T>(args)...);
  } catch (...) {  // NOLINT
  }
}

template <typename... T>
inline void raise_python_warning(PyObject* warning_type, fmt::format_string<T...> fmt,
                                 T&&... args) noexcept {
  try {
    const auto msg = fmt::format(fmt, std::forward<T>(args)...);
    [[maybe_unused]] const GilScopedAcquire gil{true};
    if (PyErr_WarnEx(warning_type, msg.c_str(), 1) < 0) {
      println_stderr_noexcept(fmt, std::forward<T>(args)...);
    }
  } catch (...) {  // NOLINT
    println_stderr_noexcept(fmt, std::forward<T>(args)...);
  }
}

template <typename... T>
inline void raise_python_user_warning(fmt::format_string<T...> fmt, T&&... args) noexcept {
  raise_python_warning(PyExc_UserWarning, fmt, std::forward<T>(args)...);
}

template <typename... T>
inline void raise_python_deprecation_warning(fmt::format_string<T...> fmt, T&&... args) noexcept {
  raise_python_warning(PyExc_DeprecationWarning, fmt, std::forward<T>(args)...);
}

template <typename... T>
inline void raise_python_runtime_warning(fmt::format_string<T...> fmt, T&&... args) noexcept {
  raise_python_warning(PyExc_RuntimeWarning, fmt, std::forward<T>(args)...);
}

template <typename T>
[[nodiscard]] inline nanobind::capsule make_capsule(std::unique_ptr<T> ptr) {
  HICTKPY_GIL_SCOPED_ACQUIRE
  nanobind::capsule owner{ptr.get(), [](void* p) noexcept {
                            delete static_cast<T*>(p);  // NOLINT
                          }};
  ptr.release();
  return owner;
}

template <typename T>
[[nodiscard]] inline nanobind::capsule make_capsule(std::unique_ptr<T> ptr, const char* name) {
  if (!name) {
    return make_capsule(std::move(ptr));
  }

  HICTKPY_GIL_SCOPED_ACQUIRE
  nanobind::capsule owner{ptr.get(), name, [](void* p) noexcept {
                            delete static_cast<T*>(p);  // NOLINT
                          }};
  ptr.release();
  return owner;
}

}  // namespace hictkpy
