// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/hic.hpp"

#include <cstdint>
#include <string>
#include <string_view>

#include "hictk/hic.hpp"
#include "hictkpy/common.hpp"

namespace py = pybind11;

namespace hictkpy::hic {
hictk::hic::File file_ctor(std::string_view path, std::int32_t resolution,
                           std::string_view matrix_type, std::string_view matrix_unit) {
  return hictk::hic::File{std::string{path}, static_cast<std::uint32_t>(resolution),
                          hictk::hic::ParseMatrixTypeStr(std::string{matrix_type}),
                          hictk::hic::ParseUnitStr(std::string{matrix_unit})};
}

py::object fetch(const hictk::hic::File &f, std::string_view range1, std::string_view range2,
                 std::string_view normalization, std::string_view count_type, bool join,
                 std::string_view query_type) {
  return file_fetch(f, range1, range2, normalization, count_type, join, query_type);
}

py::object fetch_sparse(const hictk::hic::File &f, std::string_view range1, std::string_view range2,
                        std::string_view normalization, std::string_view count_type,
                        std::string_view query_type) {
  return file_fetch_sparse(f, range1, range2, normalization, count_type, query_type);
}

py::object fetch_dense(const hictk::hic::File &f, std::string_view range1, std::string_view range2,
                       std::string_view normalization, std::string_view count_type,
                       std::string_view query_type) {
  return file_fetch_dense(f, range1, range2, normalization, count_type, query_type);
}

py::object fetch_sum(const hictk::hic::File &f, std::string_view range1, std::string_view range2,
                     std::string_view normalization, std::string_view count_type,
                     std::string_view query_type) {
  return file_fetch_sum(f, range1, range2, normalization, count_type, query_type);
}

std::int64_t fetch_nnz(const hictk::hic::File &f, std::string_view range1, std::string_view range2,
                       std::string_view query_type) {
  return file_fetch_nnz(f, range1, range2, query_type);
}
}  // namespace hictkpy::hic
