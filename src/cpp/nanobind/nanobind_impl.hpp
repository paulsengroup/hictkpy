// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Python.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <nanobind/nanobind.h>

#include <memory>
#include <utility>

#include "hictkpy/locking.hpp"

namespace hictkpy {

inline void println_stderr_noexcept(const char* msg) noexcept {
  try {
    fmt::println(stderr, FMT_STRING("{}"), msg);
  } catch (...) {  // NOLINT
  }
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
    const auto msg = fmt::format(fmt, args...);
    [[maybe_unused]] const GilScopedAcquire gil{true};
    if (PyErr_WarnEx(warning_type, msg.c_str(), 1) < 0) {
      println_stderr_noexcept(msg.c_str());
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
