// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <pybind11/pybind11.h>

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/cooler/cooler.hpp"

namespace hictkpy::cooler {

hictk::cooler::File file_ctor(std::string_view uri);

hictk::cooler::File file_ctor(std::string_view uri, const pybind11::dict &py_chroms,
                              std::uint32_t bin_size, bool overwrite_if_exists = false);
bool is_cooler(std::string_view uri);

hictk::cooler::File cooler_ctor(std::string_view uri, const pybind11::dict &py_chroms,
                                std::uint32_t bin_size, bool overwrite_if_exists = false,
                                bool float_pixels = false);
[[nodiscard]] pybind11::dict get_cooler_attrs(const hictk::cooler::File &clr);

pybind11::object fetch(const hictk::cooler::File &f, std::string_view range1,
                       std::string_view range2, std::string_view normalization,
                       std::string_view count_type, bool join, std::string_view query_type);
pybind11::object fetch_sparse(const hictk::cooler::File &f, std::string_view range1,
                              std::string_view range2, std::string_view normalization,
                              std::string_view count_type, std::string_view query_type);

pybind11::object fetch_dense(const hictk::cooler::File &f, std::string_view range1,
                             std::string_view range2, std::string_view normalization,
                             std::string_view count_type, std::string_view query_type);

pybind11::object fetch_sum(const hictk::cooler::File &f, std::string_view range1,
                           std::string_view range2, std::string_view normalization,
                           std::string_view count_type, std::string_view query_type);
std::int64_t fetch_nnz(const hictk::cooler::File &f, std::string_view range1,
                       std::string_view range2, std::string_view query_type);

}  // namespace hictkpy::cooler
