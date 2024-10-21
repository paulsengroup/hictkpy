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
#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
HICTKPY_DISABLE_WARNING_POP
// clang-format on

#include <cassert>
#include <mutex>
#include <string>

inline nanobind::module_ import_module_checked(const std::string& module_name) {
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
inline nanobind::module_ import_pyarrow_checked(int min_version_major = 16,
                                                int min_version_minor = 0,
                                                int min_version_patch = 0) {
  assert(min_version_major >= 0);
  assert(min_version_minor >= 0);
  assert(min_version_patch >= 0);

  static bool version_ok{false};
  static std::mutex mtx{};

  auto pa = import_module_checked("pyarrow");

  [[maybe_unused]] const auto lck = std::scoped_lock(mtx);
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
      // Poor man's formatting
      error_msg = "pyarrow ";
      for (const auto& tok : version) {
        error_msg += tok + ".";
      }
      if (version.size() > 1) {
        error_msg.pop_back();
      }
      error_msg +=
          " is too old to be used with hictkpy: please "
          "install a compatible version with: pip install 'hictkpy[pyarrow]'";
      throw nanobind::import_error(error_msg.c_str());
    }
  } catch (const nanobind::builtin_exception&) {
    throw;
  } catch (const std::exception& e) {
    error_msg = "unable to parse pyarrow version: ";
    error_msg += e.what();
    error_msg +=
        ". Assuming pyarrow's version is not compatible: please "
        "install a compatible version of pyarrow with: pip install 'hictkpy[pyarrow]'";
    throw nanobind::import_error(error_msg.c_str());
  } catch (...) {
    throw nanobind::import_error(
        "unable to parse pyarrow version: Assuming pyarrow's version is not compatible: please "
        "install a compatible version of pyarrow with: pip install 'hictkpy[pyarrow]'");
  }

  assert(version_ok);
  return pa;
}
