// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Python.h>

#include <memory>

namespace arrow {
class ChunkedArray;
class Table;
class Schema;

}  // namespace arrow

namespace hictkpy {

[[nodiscard]] PyObject* ExportArrowSchemaPyCapsule(const arrow::Schema& schema_in);
[[nodiscard]] PyObject* ExportArrowArrayStreamPyCapsule(
    std::shared_ptr<arrow::ChunkedArray> column);

[[nodiscard]] PyObject* export_pyarrow_table(std::shared_ptr<arrow::Table>& arrow_table);
[[nodiscard]] PyObject* export_pandas_table(PyObject* pyarrow_table);

}  // namespace hictkpy
