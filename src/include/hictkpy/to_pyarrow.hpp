// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include "hictkpy/nanobind.hpp"

namespace arrow {
class Table;

}  // namespace arrow

namespace hictkpy {

[[nodiscard]] nanobind::object export_pyarrow_table(std::shared_ptr<arrow::Table> arrow_table);

}  // namespace hictkpy
