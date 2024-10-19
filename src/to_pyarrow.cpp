// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/to_pyarrow.hpp"

#include <Python.h>
#include <arrow/c/bridge.h>
#include <arrow/python/pyarrow.h>
#include <arrow/table.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <stdexcept>

namespace hictkpy {

// NOLINTBEGIN(*-owning-memory,*-no-malloc)
static void ReleaseArrowSchemaPyCapsule(PyObject* capsule) {
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#pycapsule-standard
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html
  auto* schema = static_cast<ArrowSchema*>(PyCapsule_GetPointer(capsule, "arrow_schema"));
  if (schema->release) {
    schema->release(schema);
  }
  free(schema);
}

static void ReleaseArrowArrayStreamPyCapsule(PyObject* capsule) {
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#pycapsule-standard
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html
  auto* stream =
      static_cast<ArrowArrayStream*>(PyCapsule_GetPointer(capsule, "arrow_array_stream"));
  if (stream->release) {
    stream->release(stream);
  }
  free(stream);
}

PyObject* ExportArrowSchemaPyCapsule(const arrow::Schema& schema_in) {
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#pycapsule-standard
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html
  auto* schema_out = static_cast<ArrowSchema*>(malloc(sizeof(ArrowSchema)));
  if (!schema_out) {
    throw std::bad_alloc();
  }

  const auto status = arrow::ExportSchema(schema_in, schema_out);
  if (!status.ok()) {
    free(schema_out);
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to export arrow::Schema as ArrowSchema: {}"), status.message()));
  }

  return PyCapsule_New(schema_out, "arrow_schema", ReleaseArrowSchemaPyCapsule);
}

PyObject* ExportArrowArrayStreamPyCapsule(std::shared_ptr<arrow::ChunkedArray> column) {
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#pycapsule-standard
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html
  auto* array_stream = static_cast<ArrowArrayStream*>(malloc(sizeof(ArrowArrayStream)));
  if (!array_stream) {
    throw std::bad_alloc();
  }

  const auto status = arrow::ExportChunkedArray(std::move(column), array_stream);
  if (!status.ok()) {
    free(array_stream);
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to export arrow::ChunkedArray as ArrowArrayStream: {}"),
                    status.message()));
  }

  return PyCapsule_New(array_stream, "arrow_array_stream", ReleaseArrowArrayStreamPyCapsule);
}

PyObject* export_pyarrow_table(std::shared_ptr<arrow::Table>& arrow_table) {
  if (!arrow_table) {
    return nullptr;
  }

  std::vector<PyObject*> columns_py(static_cast<std::size_t>(arrow_table->num_columns()));
  auto* schema = ExportArrowSchemaPyCapsule(*arrow_table->schema());

  std::vector<std::shared_ptr<arrow::ChunkedArray>> columns(columns_py.size());
  std::copy(arrow_table->columns().begin(), arrow_table->columns().end(), columns.begin());
  arrow_table.reset();

  for (std::size_t i = 0; i < columns.size(); ++i) {
    // NOLINTNEXTLINE(*-pro-bounds-pointer-arithmetic)
    columns_py[i] = ExportArrowArrayStreamPyCapsule(std::move(columns[i]));
  }

  return table;
}
// NOLINTEND(*-owning-memory,*-no-malloc)

}  // namespace hictkpy
