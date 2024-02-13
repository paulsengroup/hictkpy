// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <nanobind/nanobind.h>

#include <string>
#include <string_view>

#include "hictkpy/multires_file.hpp"

namespace nb = nanobind;

namespace hictkpy::multires_file {

void ctor(hictk::cooler::MultiResFile* fp, std::string_view path) {
  new (fp) hictk::cooler::MultiResFile(std::string{path});
}

std::string repr(const hictk::cooler::MultiResFile& mclr) {
  return fmt::format(FMT_STRING("MultiResFile({})"), mclr.path());
}

nb::dict get_attrs(const hictk::cooler::MultiResFile& mclr) {
  nb::dict py_attrs;

  py_attrs["format"] = mclr.attributes().format;
  py_attrs["format_version"] = mclr.attributes().format_version;
  py_attrs["bin_type"] = mclr.attributes().bin_type;

  return py_attrs;
}

hictk::File getitem(const hictk::cooler::MultiResFile& mclr, std::uint32_t resolution) {
  return hictk::File(mclr.open(resolution));
}

}  // namespace hictkpy::multires_file
