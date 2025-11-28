// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/bin_table.hpp"

#include <arrow/array.h>
#include <arrow/buffer.h>
#include <arrow/builder.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <hictk/bin.hpp>
#include <hictk/bin_table.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/reference.hpp>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/reference.hpp"
#include "hictkpy/table.hpp"

namespace nb = nanobind;

namespace hictkpy {

BinTable::BinTable(std::shared_ptr<const hictk::BinTable> bins) noexcept : _bins(std::move(bins)) {
  if (!_bins) {
    _bins = std::make_shared<const hictk::BinTable>();
  }
}

BinTable::BinTable(hictk::BinTable bins)
    : BinTable(std::make_shared<const hictk::BinTable>(std::move(bins))) {}

BinTable::BinTable(const ChromosomeDict& chromosomes, std::uint32_t resolution)
    : BinTable(hictk::BinTable{chromosome_dict_to_reference(chromosomes), resolution}) {}

[[noreturn]] static void throw_except_failed_to_parse_bins_df(const std::exception& e) {
  throw std::runtime_error(fmt::format(
      FMT_STRING(
          "Unable to fetch bins from the given object. Please make sure the given object is a "
          "pandas.DataFrame with columns [\"chrom\", \"start\", \"end\"]. Underlying error: {}"),
      e.what()));
}

[[noreturn]] static void throw_except_failed_to_parse_bins_df() {
  throw_except_failed_to_parse_bins_df(std::runtime_error{"unknown error"});
}

[[nodiscard]] static hictk::Reference get_reference_from_bins_df(const nb::object& df) {
  try {
    return chromosome_dict_to_reference(
        nb::cast<ChromosomeDict>(df.attr("groupby")("chrom", nb::arg("observed") = true)
                                     .attr("__getitem__")("end")
                                     .attr("max")()
                                     .attr("to_dict")()));
  } catch (const std::exception& e) {
    throw_except_failed_to_parse_bins_df(e);
  } catch (...) {
    throw_except_failed_to_parse_bins_df();
  }
}

template <typename I>
[[nodiscard]] static std::vector<I> get_std_vect_from_bins_df(const nb::object& df,
                                                              std::string_view col_name) {
  static_assert(std::is_arithmetic_v<I>);
  try {
    return nb::cast<std::vector<I>>(df.attr("__getitem__")(col_name));
  } catch (const std::exception& e) {
    throw_except_failed_to_parse_bins_df(e);
  } catch (...) {
    throw_except_failed_to_parse_bins_df();
  }
}

BinTable::BinTable(const nanobind::object& df)
    : BinTable(get_reference_from_bins_df(df),
               get_std_vect_from_bins_df<std::uint32_t>(df, "start"),
               get_std_vect_from_bins_df<std::uint32_t>(df, "end")) {}

BinTable::BinTable(hictk::Reference chroms, const std::vector<std::uint32_t>& start_pos,
                   const std::vector<std::uint32_t>& end_pos)
    : BinTable(hictk::BinTable{std::move(chroms), start_pos, end_pos}) {}

const hictk::Reference& BinTable::chromosomes() const noexcept { return _bins->chromosomes(); }

std::uint32_t BinTable::resolution() const noexcept { return _bins->resolution(); }

std::string_view BinTable::type() const noexcept {
  return _bins->type() == hictk::BinTable::Type::fixed ? "fixed" : "variable";
}

std::size_t BinTable::size() const noexcept { return _bins->size(); }

std::string BinTable::repr() const {
  return fmt::format(FMT_STRING("BinTable(num_chroms={}; bin_size={};)"),
                     _bins->chromosomes().size(),
                     _bins->resolution() == 0 ? "variable" : fmt::to_string(_bins->resolution()));
}

[[nodiscard]] static std::shared_ptr<arrow::DataType> chrom_dict() {
  return dictionary(arrow::int32(), arrow::utf8());
}

[[nodiscard]] static std::shared_ptr<arrow::Schema> make_bin_table_schema(
    bool include_bin_id = false) {
  arrow::FieldVector fields{};
  fields.reserve(include_bin_id ? 4 : 3);

  if (include_bin_id) {
    fields.emplace_back(arrow::field("bin_id", arrow::uint64(), false));
  }

  fields.emplace_back(arrow::field("chrom", chrom_dict(), false));
  fields.emplace_back(arrow::field("start", arrow::uint32(), false));
  fields.emplace_back(arrow::field("end", arrow::uint32(), false));

  return arrow::schema(fields);
}

[[nodiscard]] static std::shared_ptr<arrow::Array> chrom_names_to_arrow(
    const std::vector<std::string>& chrom_names) {
  arrow::StringBuilder chrom_name_builder{};
  auto status = chrom_name_builder.AppendValues(chrom_names);
  if (!status.ok()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to construct a table of pixels: {}"), status.message()));
  }

  auto result = chrom_name_builder.Finish();
  if (!result.ok()) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to construct a table of pixels: {}"),
                                         result.status().message()));
  }

  return result.MoveValueUnsafe();
}

[[nodiscard]] static std::shared_ptr<arrow::Table> bin_table_to_arrow(
    const std::vector<std::string>& chrom_names, std::vector<std::int32_t> chrom_ids,
    std::vector<std::uint32_t> start_pos, std::vector<std::uint32_t> end_pos,
    std::optional<std::vector<std::uint64_t>> bin_ids) {
  const auto num_bins = static_cast<std::int64_t>(chrom_ids.size());

  assert(chrom_ids.size() == start_pos.size());
  assert(chrom_ids.size() == end_pos.size());
  if (bin_ids.has_value()) {
    assert(chrom_ids.size() == bin_ids->size());
  }

  if (num_bins == 0) {
    return arrow::Table::Make(make_bin_table_schema(bin_ids.has_value()),
                              std::vector<std::shared_ptr<arrow::Array>>{}, 0);
  }

  assert(!chrom_names.empty());

  std::shared_ptr<arrow::UInt64Array> bin_ids_data{};
  if (bin_ids.has_value()) {
    bin_ids_data = std::make_shared<arrow::UInt64Array>(
        num_bins, arrow::Buffer::FromVector(std::move(*bin_ids)), nullptr, 0, 0);
  }

  auto chrom_data = std::make_shared<arrow::DictionaryArray>(
      chrom_dict(),
      std::make_shared<arrow::Int32Array>(num_bins, arrow::Buffer::FromVector(std::move(chrom_ids)),
                                          nullptr, 0, 0),
      chrom_names_to_arrow(chrom_names));

  auto start_data = std::make_shared<arrow::UInt32Array>(
      num_bins, arrow::Buffer::FromVector(std::move(start_pos)), nullptr, 0, 0);
  auto end_data = std::make_shared<arrow::UInt32Array>(
      num_bins, arrow::Buffer::FromVector(std::move(end_pos)), nullptr, 0, 0);

  auto schema = make_bin_table_schema(!!bin_ids_data);

  std::vector<std::shared_ptr<arrow::Array>> data{};
  data.reserve(!!bin_ids_data ? 4 : 3);

  if (bin_ids_data) {
    data.emplace_back(std::move(bin_ids_data));
  }

  data.emplace_back(std::move(chrom_data));
  data.emplace_back(std::move(start_data));
  data.emplace_back(std::move(end_data));

  return arrow::Table::Make(std::move(schema), data);
}

[[nodiscard]] static nb::object make_bin_table_pyarrow(
    const std::vector<std::string>& chrom_names, std::vector<std::int32_t> chrom_ids,
    std::vector<std::uint32_t> start_pos, std::vector<std::uint32_t> end_pos,
    std::optional<std::vector<std::uint64_t>> bin_ids = {}) {
  return export_pyarrow_table(bin_table_to_arrow(chrom_names, std::move(chrom_ids),
                                                 std::move(start_pos), std::move(end_pos),
                                                 std::move(bin_ids)));
}

[[nodiscard]] static nb::object make_bin_table_df(
    const std::vector<std::string>& chrom_names, std::vector<std::int32_t> chrom_ids,
    std::vector<std::uint32_t> start_pos, std::vector<std::uint32_t> end_pos,
    std::optional<std::vector<std::uint64_t>> bin_ids = {}) {
  auto table = make_bin_table_pyarrow(chrom_names, std::move(chrom_ids), std::move(start_pos),
                                      std::move(end_pos), std::move(bin_ids));

  HICTKPY_GIL_SCOPED_ACQUIRE
  return nb::cast(table.attr("to_pandas")(nb::arg("self_destruct") = true),
                  nb::rv_policy::take_ownership);
}

nb::object BinTable::bin_ids_to_coords(std::vector<std::uint64_t> bin_ids) const {
  std::vector<std::int32_t> chrom_ids(bin_ids.size());
  std::vector<std::uint32_t> start_pos(bin_ids.size());
  std::vector<std::uint32_t> end_pos(bin_ids.size());

  std::visit(
      [&](const auto& bins) {
        for (std::size_t i = 0; i < bin_ids.size(); ++i) {
          const auto bin = bins.at(bin_ids[i]);
          chrom_ids[i] = static_cast<std::int32_t>(bin.chrom().id());
          start_pos[i] = bin.start();
          end_pos[i] = bin.end();
        }
      },
      _bins->get());

  return make_bin_table_df(chrom_names(), std::move(chrom_ids), std::move(start_pos),
                           std::move(end_pos), std::move(bin_ids));
}

hictk::Bin BinTable::bin_id_to_coord(std::uint64_t bin_id) const { return _bins->at(bin_id); }

hictk::Bin BinTable::coord_to_bin(std::string_view chrom, std::uint32_t pos) const {
  return _bins->at(chrom, pos);
}

nanobind::object BinTable::coords_to_bins(const std::vector<std::string>& chroms,
                                          const std::vector<std::uint32_t>& positions) const {
  if (chroms.size() != positions.size()) {
    throw std::runtime_error("chroms and positions should have the same size");
  }
  std::vector<std::uint64_t> bin_ids(chroms.size());
  std::vector<std::int32_t> chrom_ids(chroms.size());
  std::vector<std::uint32_t> start_pos(chroms.size());
  std::vector<std::uint32_t> end_pos(chroms.size());

  std::visit(
      [&](const auto& bins) {
        for (std::size_t i = 0; i < chroms.size(); ++i) {
          const auto bin = bins.at(chroms[i], positions[i]);
          bin_ids[i] = bin.id();
          chrom_ids[i] = static_cast<std::int32_t>(bin.chrom().id());
          start_pos[i] = bin.start();
          end_pos[i] = bin.end();
        }
      },
      _bins->get());

  return make_bin_table_df(chrom_names(), std::move(chrom_ids), std::move(start_pos),
                           std::move(end_pos), std::move(bin_ids));
}

std::int64_t BinTable::coord_to_bin_id(std::string_view chrom, std::uint32_t pos) const {
  return static_cast<std::int64_t>(_bins->at(chrom, pos).id());
}

auto BinTable::coords_to_bin_ids(const std::vector<std::string>& chroms,
                                 const std::vector<std::uint32_t>& positions) const -> BinIDsVec {
  auto np = import_module_checked("numpy");

  if (chroms.size() != positions.size()) {
    throw std::runtime_error("chroms and positions should have the same size");
  }

  auto np_array =
      np.attr("empty")(static_cast<std::int64_t>(chroms.size()), nb::arg("dtype") = "int");
  auto buffer = nb::cast<BinIDsVec>(np_array);
  auto bin_ids = buffer.view();

  std::visit(
      [&](const auto& bins) {
        for (std::size_t i = 0; i < chroms.size(); ++i) {
          bin_ids(static_cast<std::int64_t>(i)) =
              static_cast<std::int64_t>(bins.at(chroms[i], positions[i]).id());
        }
      },
      _bins->get());

  return buffer;
}

[[nodiscard]] static nb::object make_bg2_pixels_df(const std::vector<std::string>& chrom_names,
                                                   std::vector<std::int32_t> chrom1_ids,
                                                   std::vector<std::uint32_t> start1_pos,
                                                   std::vector<std::uint32_t> end1_pos,
                                                   std::vector<std::int32_t> chrom2_ids,
                                                   std::vector<std::uint32_t> start2_pos,
                                                   std::vector<std::uint32_t> end2_pos) {
  auto pd = import_module_checked("pandas");

  nb::list dfs{};
  dfs.append(make_bin_table_df(chrom_names, std::move(chrom1_ids), std::move(start1_pos),
                               std::move(end1_pos)));
  dfs.append(make_bin_table_df(chrom_names, std::move(chrom2_ids), std::move(start2_pos),
                               std::move(end2_pos)));

  auto df = pd.attr("concat")(dfs, nb::arg("axis") = "columns", nb::arg("ignore_index") = true,
                              nb::arg("copy") = false);

  static const std::vector<std::string_view> col_names{"chrom1", "start1", "end1",
                                                       "chrom2", "start2", "end2"};
  df.attr("columns") = col_names;
  return df;
}

nb::object BinTable::merge_coords(nb::object df) const {
  check_pyarrow_is_importable();
  auto pd = import_module_checked("pandas");

  // TODO avoid potentially copying data here
  using Buffer64T = nb::ndarray<nb::numpy, nb::shape<-1>, std::int64_t>;
  auto bin1_ids = nb::cast<Buffer64T>(df.attr("__getitem__")("bin1_id").attr("to_numpy")());
  auto bin2_ids = nb::cast<Buffer64T>(df.attr("__getitem__")("bin2_id").attr("to_numpy")());

  const auto n = bin1_ids.size();

  std::vector<std::int32_t> chrom1_ids(n);
  std::vector<std::uint32_t> starts1(n);
  std::vector<std::uint32_t> ends1(n);
  std::vector<std::int32_t> chrom2_ids(n);
  std::vector<std::uint32_t> starts2(n);
  std::vector<std::uint32_t> ends2(n);

  std::visit(
      [&](const auto& bins) {
        for (std::size_t i = 0; i < n; ++i) {
          const auto bin1 =
              bins.at(static_cast<std::uint64_t>(bin1_ids(static_cast<std::int64_t>(i))));
          chrom1_ids[i] = static_cast<std::int32_t>(bin1.chrom().id());
          starts1[i] = bin1.start();
          ends1[i] = bin1.end();

          const auto bin2 =
              bins.at(static_cast<std::uint64_t>(bin2_ids(static_cast<std::int64_t>(i))));
          chrom2_ids[i] = static_cast<std::int32_t>(bin2.chrom().id());
          starts2[i] = bin2.start();
          ends2[i] = bin2.end();
        }
      },
      _bins->get());

  auto coord_df =
      make_bg2_pixels_df(chrom_names(), std::move(chrom1_ids), std::move(starts1), std::move(ends1),
                         std::move(chrom2_ids), std::move(starts2), std::move(ends2));

  std::vector<nb::any> index_names{};
  try {
    auto name = nb::cast<std::optional<nb::any>>(df.attr("index").attr("name"));
    index_names.emplace_back(name.value_or(nb::cast<nb::any>(nb::str{"index"})));
  } catch (const std::bad_cast&) {
    index_names = nb::cast<std::vector<nb::any>>(df.attr("index").attr("name"));
  }

  auto col_names = index_names;
  for (auto&& name : nb::cast<std::vector<nb::any>>(df.attr("columns").attr("tolist")())) {
    col_names.emplace_back(name);
  }
  for (auto&& name : nb::cast<std::vector<nb::any>>(coord_df.attr("columns").attr("tolist")())) {
    col_names.emplace_back(name);
  }

  nb::list dfs{};
  dfs.append(df.attr("reset_index")());
  dfs.append(std::move(coord_df));

  df = pd.attr("concat")(dfs, nb::arg("axis") = "columns", nb::arg("ignore_index") = true,
                         nb::arg("copy") = false);
  df.attr("columns") = col_names;
  df.attr("set_index")(index_names, nb::arg("inplace") = true);
  return df;
}

nb::iterator BinTable::make_iterable() const {
  return std::visit(
      [](const auto& bins) {
        return nb::make_iterator(nb::type<hictkpy::BinTable>(), "BinTableIterator", bins.begin(),
                                 bins.end());
      },
      _bins->get());
}

static auto compute_num_bins(const hictk::BinTable& bins, const hictk::GenomicInterval& query) {
  if (!query) {
    return bins.size();
  }
  return static_cast<std::size_t>(std::visit(
      [&](const auto& bins_) {
        const auto [first_bin, last_bin] = bins_.find_overlap(query);
        return std::distance(first_bin, last_bin);
      },
      bins.get()));
}

nb::object BinTable::to_arrow(std::optional<std::string_view> range,
                              std::string_view query_type) const {
  const auto qt =
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED;
  const auto query = !range.has_value() ? hictk::GenomicInterval{}
                                        : hictk::GenomicInterval::parse(
                                              _bins->chromosomes(), std::string{range.value()}, qt);

  const auto n = compute_num_bins(*_bins, query);

  std::vector<std::uint64_t> bin_ids(n);
  std::vector<std::int32_t> chrom_ids(n);
  std::vector<std::uint32_t> starts(n);
  std::vector<std::uint32_t> ends(n);

  const auto chrom_id_offset = static_cast<std::uint32_t>(_bins->chromosomes().at(0).is_all());

  std::visit(
      [&](const auto& bins) {
        const auto [first_bin, last_bin] = !range.has_value()
                                               ? std::make_pair(bins.begin(), bins.end())
                                               : bins.find_overlap(query);
        std::size_t i = 0;
        std::for_each(first_bin, last_bin, [&](const auto& bin) {
          bin_ids[i] = bin.id();
          chrom_ids[i] = static_cast<std::int32_t>(bin.chrom().id() - chrom_id_offset);
          starts[i] = bin.start();
          ends[i] = bin.end();
          ++i;
        });
      },
      _bins->get());

  return make_bin_table_pyarrow(chrom_names(false), std::move(chrom_ids), std::move(starts),
                                std::move(ends), std::move(bin_ids));
}

nb::object BinTable::to_pandas(std::optional<std::string_view> range,
                               std::string_view query_type) const {
  auto table = std::make_optional(to_arrow(range, query_type));
  HICTKPY_GIL_SCOPED_ACQUIRE
  auto df = nb::cast(table->attr("to_pandas")(nb::arg("self_destruct") = true),
                     nb::rv_policy::take_ownership);
  table.reset();
  return df;
}

std::shared_ptr<const hictk::BinTable> BinTable::get() const noexcept { return _bins; }

std::vector<std::string> BinTable::chrom_names(bool include_ALL) const {
  std::vector<std::string> chrom_names;
  chrom_names.reserve(_bins->chromosomes().size());
  for (const auto& chrom : _bins->chromosomes()) {
    if (!include_ALL && chrom.is_all()) {
      continue;
    }

    chrom_names.emplace_back(chrom.name());
  }

  return chrom_names;
}

static void declare_bin_class(nb::module_& m) {
  nb::class_<hictk::Bin>(m, "Bin", "Class representing a genomic Bin (i.e., a BED interval).")
      .def_prop_ro(
          "id", [](const hictk::Bin& b) { return b.id(); }, "Get the bin ID.")
      .def_prop_ro(
          "rel_id", [](const hictk::Bin& b) { return b.rel_id(); },
          "Get the relative bin ID "
          "(i.e., the ID that uniquely identifies a bin within a chromosome).")
      .def_prop_ro(
          "chrom", [](const hictk::Bin& b) { return b.chrom().name(); },
          "Get the name of the chromosome to which the Bin refers to.")
      .def_prop_ro(
          "start", [](const hictk::Bin& b) { return b.start(); }, "Get the Bin start position.")
      .def_prop_ro(
          "end", [](const hictk::Bin& b) { return b.end(); }, "Get the Bin end position.")
      .def("__repr__",
           [](const hictk::Bin& b) {
             return fmt::format(FMT_COMPILE("id={}; rel_id={}; chrom={}; start={}; end={}"), b.id(),
                                b.rel_id(), b.chrom().name(), b.start(), b.end());
           })
      .def("__str__", [](const hictk::Bin& b) {
        return fmt::format(FMT_COMPILE("{}\t{}\t{}"), b.chrom().name(), b.start(), b.end());
      });
}

void BinTable::bind(nb::module_& m) {
  declare_bin_class(m);

  auto bt = nb::class_<BinTable>(m, "BinTable", "Class representing a table of genomic bins.");

  bt.def(nb::init<ChromosomeDict, std::uint32_t>(), nb::arg("chroms"), nb::arg("resolution"),
         "Construct a table of bins given a dictionary mapping chromosomes to their sizes and a "
         "resolution.");

  bt.def(nb::init<nb::object>(), nb::arg("bins"),
         "Construct a table of bins from a pandas.DataFrame with columns [\"chrom\", \"start\", "
         "\"end\"].",
         nb::sig("def __init__(self, bins: pandas.DataFrame) -> None"));

  bt.def("__repr__", &BinTable::repr, nb::rv_policy::move);

  bt.def("chromosomes", &get_chromosomes_from_object<hictkpy::BinTable>,
         nb::arg("include_ALL") = false,
         "Get the chromosome sizes as a dictionary mapping names to sizes.",
         nb::rv_policy::take_ownership);

  bt.def("resolution", &BinTable::resolution,
         "Get the bin size for the bin table. "
         "Return 0 in case the bin table has a variable bin size.");

  bt.def("type", &BinTable::type,
         "Get the type of table underlying the BinTable object (i.e. fixed or variable).",
         nb::rv_policy::move);

  bt.def("__len__", &BinTable::size, "Get the number of bins in the bin table.");

  bt.def("__iter__", &BinTable::make_iterable, nb::keep_alive<0, 1>(),
         nb::sig("def __iter__(self) -> hictkpy.BinTableIterator"),
         "Implement iter(self). The resulting iterator yields objects of type hictkpy.Bin.",
         nb::rv_policy::take_ownership);

  bt.def("get", &BinTable::bin_id_to_coord, nb::arg("bin_id"),
         "Get the genomic coordinate given a bin ID.",
         nb::sig("def get(self, bin_id: int) -> hictkpy.Bin"), nb::rv_policy::move);
  bt.def("get", &BinTable::bin_ids_to_coords, nb::arg("bin_ids"),
         "Get the genomic coordinates given a sequence of bin IDs. "
         "Genomic coordinates are returned as a pandas.DataFrame with columns [\"chrom\", "
         "\"start\", \"end\"].",
         nb::sig("def get(self, bin_ids: collections.abc.Sequence[int]) -> pandas.DataFrame"),
         nb::rv_policy::take_ownership);

  bt.def("get", &BinTable::coord_to_bin, nb::arg("chrom"), nb::arg("pos"),
         "Get the bin overlapping the given genomic coordinate.",
         nb::sig("def get(self, chrom: str, pos: int) -> hictkpy.Bin"), nb::rv_policy::move);
  bt.def("get", &BinTable::coords_to_bins, nb::arg("chroms"), nb::arg("pos"),
         "Get the bins overlapping the given genomic coordinates. "
         "Bins are returned as a pandas.DataFrame with columns [\"chrom\", "
         "\"start\", \"end\"].",
         nb::sig("def get(self, chroms: collections.abc.Sequence[str], pos: "
                 "collections.abc.Sequence[int]) -> pandas.DataFrame"),
         nb::rv_policy::take_ownership);

  bt.def("get_id", &BinTable::coord_to_bin_id, nb::arg("chrom"), nb::arg("pos"),
         "Get the ID of the bin overlapping the given genomic coordinate.");
  bt.def("get_ids", &BinTable::coords_to_bin_ids, nb::arg("chroms"), nb::arg("pos"),
         "Get the IDs of the bins overlapping the given genomic coordinates.",
         nb::rv_policy::take_ownership);

  bt.def("merge", &BinTable::merge_coords, nb::arg("df"),
         "Merge genomic coordinates corresponding to the given bin identifiers. "
         "Bin identifiers should be provided as a pandas.DataFrame with columns \"bin1_id\" and "
         "\"bin2_id\". "
         "Genomic coordinates are returned as a pandas.DataFrame containing the same data as the "
         "DataFrame given as input, plus columns [\"chrom1\", \"start1\", \"end1\", \"chrom2\", "
         "\"start2\", \"end2\"].",
         nb::sig("def merge(self, df: pandas.DataFrame) -> pandas.DataFrame"),
         nb::rv_policy::take_ownership);

  bt.def("to_arrow", &BinTable::to_arrow, nb::call_guard<nb::gil_scoped_release>(),
         nb::arg("range") = nb::none(), nb::arg("query_type") = "UCSC",
         "Return the bins in the BinTable as a pyarrow.Table. The optional \"range\" parameter "
         "can be used to only fetch a subset of the bins in the BinTable.",
         nb::sig("def to_arrow(self, range: str | None = None, query_type: str = 'UCSC') -> "
                 "pyarrow.Table"),
         nb::rv_policy::take_ownership);

  bt.def("to_pandas", &BinTable::to_pandas, nb::call_guard<nb::gil_scoped_release>(),
         nb::arg("range") = nb::none(), nb::arg("query_type") = "UCSC",
         "Return the bins in the BinTable as a pandas.DataFrame. The optional \"range\" parameter "
         "can be used to only fetch a subset of the bins in the BinTable.",
         nb::sig("def to_pandas(self, range: str | None = None, query_type: str = 'UCSC') -> "
                 "pandas.DataFrame"),
         nb::rv_policy::take_ownership);

  bt.def("to_df", &BinTable::to_pandas, nb::call_guard<nb::gil_scoped_release>(),
         nb::arg("range") = nb::none(), nb::arg("query_type") = "UCSC", "Alias to to_pandas().",
         nb::sig("def to_df(self, range: str | None = None, query_type: str = 'UCSC') -> "
                 "pandas.DataFrame"),
         nb::rv_policy::take_ownership);
}

}  // namespace hictkpy
