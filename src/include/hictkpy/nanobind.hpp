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

#include <string>

inline nanobind::module_ import_module_checked(const std::string& module_name) {
  try {
    return nanobind::module_::import_(module_name.c_str());
  } catch (nanobind::python_error& e) {
    // NOLINTNEXTLINE(*-pro-type-vararg)
    nanobind::raise_from(e, PyExc_ModuleNotFoundError,
                         "To enable %s support, please install %s with: pip install 'hictkpy[%s]'",
                         module_name.c_str(), module_name.c_str(), module_name.c_str());
  }
}
