// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/to_pyarrow.hpp"

#include <Python.h>
#include <arrow/c/abi.h>
#include <arrow/c/bridge.h>
#include <arrow/table.h>
#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "hictkpy/nanobind.hpp"

namespace nb = nanobind;

namespace hictkpy {

// https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#pycapsule-standard
// https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html

// NOLINTBEGIN(*-owning-memory,*-no-malloc)
static void release_arrow_schemaPyCapsule(PyObject* capsule) {
  auto* schema = static_cast<ArrowSchema*>(PyCapsule_GetPointer(capsule, "arrow_schema"));
  if (schema->release) {
    schema->release(schema);
  }
  free(schema);
}

static void release_arrow_array_streamPyCapsule(PyObject* capsule) {
  auto* stream =
      static_cast<ArrowArrayStream*>(PyCapsule_GetPointer(capsule, "arrow_array_stream"));
  if (stream->release) {
    stream->release(stream);
  }
  free(stream);
}

[[nodiscard]] static nb::object export_arrow_schema(const arrow::Schema& schema_in) {
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#arrowschema-export
  auto* schema = static_cast<ArrowSchema*>(malloc(sizeof(ArrowSchema)));
  if (!schema) {
    throw std::bad_alloc();
  }

  const auto status = arrow::ExportSchema(schema_in, schema);
  if (!status.ok()) {
    free(schema);
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to export arrow::Schema as ArrowSchema: {}"), status.message()));
  }

  auto* capsule =
      PyCapsule_New(static_cast<void*>(schema), "arrow_schema", release_arrow_schemaPyCapsule);
  if (!capsule) {
    if (schema->release) {
      schema->release(schema);
    }
    free(schema);
    throw std::runtime_error("Failed to create PyCapsule for arrow_schema");
  }

  auto obj = nb::module_::import_("types").attr("SimpleNamespace")();
  obj.attr("__setattr__")("__arrow_c_schema__",
                          nb::cpp_function([ptr = capsule]() { return nb::handle{ptr}; }));
  return obj;
}

[[nodiscard]] static nb::object export_arrow_array_stream(
    std::shared_ptr<arrow::ChunkedArray> column) {
  // https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#arrowstream-export
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

  auto* capsule = PyCapsule_New(static_cast<void*>(array_stream), "arrow_array_stream",
                                release_arrow_array_streamPyCapsule);
  if (!capsule) {
    if (array_stream->release) {
      array_stream->release(array_stream);
    }
    free(array_stream);
    throw std::runtime_error("Failed to create PyCapsule for arrow_array_stream");
  }

  auto obj = nb::module_::import_("types").attr("SimpleNamespace")();
  obj.attr("__setattr__")(
      "__arrow_c_stream__",
      nb::cpp_function(
          [ptr = capsule]([[maybe_unused]] const nb::any& _) { return nb::handle{ptr}; },
          nb::arg("requested_schema") = nb::none()));
  return obj;
}

nb::object export_pyarrow_table(std::shared_ptr<arrow::Table> arrow_table) {
  assert(arrow_table);

  const auto pa = import_pyarrow_checked();

  std::vector<nb::object> columns_py(static_cast<std::size_t>(arrow_table->num_columns()));
  std::vector<std::shared_ptr<arrow::ChunkedArray>> columns(columns_py.size());
  std::copy(arrow_table->columns().begin(), arrow_table->columns().end(), columns.begin());

  auto schema = pa.attr("schema")(export_arrow_schema(*arrow_table->schema()));
  arrow_table.reset();

  for (std::size_t i = 0; i < columns.size(); ++i) {
    // NOLINTNEXTLINE(*-pro-bounds-pointer-arithmetic)
    columns_py[i] = pa.attr("chunked_array")(export_arrow_array_stream(std::move(columns[i])));
  }

  return nb::cast(pa.attr("Table").attr("from_arrays")(columns_py, nb::arg("schema") = schema),
                  nb::rv_policy::take_ownership);
}
// NOLINTEND(*-owning-memory,*-no-malloc)

}  // namespace hictkpy
