// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/type_fwd.h>

#include <memory>

#include "hictkpy/nanobind.hpp"
#include "hictkpy/table.hpp"

namespace hictkpy::internal {

[[nodiscard]] NumericDtype infer_count_type(const std::shared_ptr<arrow::Table>& df);

[[noreturn]] void raise_invalid_dict_format();
[[noreturn]] void raise_invalid_table_format();

[[nodiscard]] PyArrowTable make_table(const nanobind::dict& columns);

}  // namespace hictkpy::internal
