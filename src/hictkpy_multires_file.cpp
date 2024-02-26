// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <string>
#include <string_view>

#include "hictkpy/multires_file.hpp"


namespace hictkpy::multires_file {

void ctor(hictk::MultiResFile* fp, std::string_view path) {
  new (fp) hictk::MultiResFile(std::string{path});
}

std::string repr(const hictk::MultiResFile& mrf) {
  return fmt::format(FMT_STRING("MultiResFile({})"), mrf.path());
}

}  // namespace hictkpy::multires_file
