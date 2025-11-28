// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Python.h>
#include <fmt/format.h>
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

#include <memory>
#include <string>

namespace hictkpy {

[[nodiscard]] nanobind::module_ import_module_checked(const std::string& module_name);

// NOLINTNEXTLINE(*-avoid-magic-numbers)
[[nodiscard]] nanobind::module_ import_pyarrow_checked(int min_version_major = 16,
                                                       int min_version_minor = 0,
                                                       int min_version_patch = 0);

void check_module_is_importable(const std::string& module_name);
// NOLINTNEXTLINE(*-avoid-magic-numbers)
void check_pyarrow_is_importable(int min_version_major = 16, int min_version_minor = 0,
                                 int min_version_patch = 0);

void println_stderr_noexcept(const char* msg) noexcept;

template <typename... T>
void println_stderr_noexcept(fmt::format_string<T...> fmt, T&&... args) noexcept;

template <typename... T>
void raise_python_warning(PyObject* warning_type, fmt::format_string<T...> fmt,
                          T&&... args) noexcept;
template <typename... T>
void raise_python_user_warning(fmt::format_string<T...> fmt, T&&... args) noexcept;

template <typename... T>
void raise_python_deprecation_warning(fmt::format_string<T...> fmt, T&&... args) noexcept;

template <typename... T>
void raise_python_runtime_warning(fmt::format_string<T...> fmt, T&&... args) noexcept;

template <typename T>
[[nodiscard]] nanobind::capsule make_capsule(std::unique_ptr<T> ptr);

template <typename T>
[[nodiscard]] nanobind::capsule make_capsule(std::unique_ptr<T> ptr, const char* name);

}  // namespace hictkpy

#include "../../nanobind_impl.hpp"
