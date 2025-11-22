// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "hictkpy/nanobind.hpp"

namespace arrow {
class Table;

}  // namespace arrow

namespace hictkpy {

[[nodiscard]] nanobind::object export_pyarrow_table(std::shared_ptr<arrow::Table> arrow_table);

void test_import_table(const nanobind::object& df, const std::vector<std::string>& column_names);

}  // namespace hictkpy
