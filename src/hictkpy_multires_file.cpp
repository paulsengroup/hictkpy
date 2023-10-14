// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <pybind11/pybind11.h>

#include <string>
#include <string_view>

#include "hictkpy/multires_file.hpp"

namespace py = pybind11;

namespace hictkpy::multires_file {

hictk::cooler::MultiResFile ctor(std::string_view path) {
  return hictk::cooler::MultiResFile(std::string{path});
}

std::string repr(const hictk::cooler::MultiResFile& mclr) {
  return fmt::format(FMT_STRING("MultiResFile({})"), mclr.path());
}

py::dict get_attrs(const hictk::cooler::MultiResFile& mclr) {
  py::dict py_attrs;

  py_attrs["format"] = mclr.attributes().format;
  py_attrs["format_version"] = mclr.attributes().format_version;
  py_attrs["bin_type"] = mclr.attributes().bin_type;

  return py_attrs;
}

hictk::File getitem(const hictk::cooler::MultiResFile& mclr, std::uint32_t resolution) {
  return hictk::File(mclr.open(resolution));
}

}  // namespace hictkpy::multires_file
