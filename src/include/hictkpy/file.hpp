// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <string_view>

#include "hictk/file.hpp"
#include "hictkpy/common.hpp"

namespace hictkpy::file {
[[nodiscard]] hictk::File ctor(std::string_view path, std::int32_t resolution,
                               std::string_view matrix_type, std::string_view matrix_unit);

pybind11::object fetch(const hictk::File &f, std::string_view range1, std::string_view range2,
                       std::string_view normalization, std::string_view count_type, bool join,
                       std::string_view query_type);

pybind11::object fetch_sparse(const hictk::File &f, std::string_view range1,
                              std::string_view range2, std::string_view normalization,
                              std::string_view count_type, std::string_view query_type);

pybind11::object fetch_dense(const hictk::File &f, std::string_view range1, std::string_view range2,
                             std::string_view normalization, std::string_view count_type,
                             std::string_view query_type);

pybind11::object fetch_sum(const hictk::File &f, std::string_view range1, std::string_view range2,
                           std::string_view normalization, std::string_view count_type,
                           std::string_view query_type);

std::int64_t fetch_nnz(const hictk::File &f, std::string_view range1, std::string_view range2,
                       std::string_view query_type);

}  // namespace hictkpy::file
