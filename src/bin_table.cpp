// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/bin_table.hpp"

#include <arrow/array.h>
#include <arrow/builder.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <hictk/bin.hpp>
#include <hictk/bin_table.hpp>
#include <hictk/reference.hpp>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "hictkpy/reference.hpp"
#include "hictkpy/to_pyarrow.hpp"

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

[[nodiscard]] static std::shared_ptr<arrow::Schema> make_bin_table_schema() {
  arrow::FieldVector fields{};
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

[[nodiscard]] static nb::object make_bin_table_df(const std::vector<std::string>& chrom_names,
                                                  std::vector<std::int32_t> chrom_ids,
                                                  std::vector<std::uint32_t> start_pos,
                                                  std::vector<std::uint32_t> end_pos,
                                                  std::vector<std::uint64_t> bin_ids = {}) {
  const auto num_bins = static_cast<std::int64_t>(chrom_ids.size());

  auto schema = make_bin_table_schema();

  auto chrom_data = std::make_shared<arrow::DictionaryArray>(
      chrom_dict(),
      std::make_shared<arrow::Int32Array>(num_bins, arrow::Buffer::FromVector(std::move(chrom_ids)),
                                          nullptr, 0, 0),
      chrom_names_to_arrow(chrom_names));

  auto start_data = std::make_shared<arrow::UInt32Array>(
      num_bins, arrow::Buffer::FromVector(std::move(start_pos)), nullptr, 0, 0);
  auto end_data = std::make_shared<arrow::UInt32Array>(
      num_bins, arrow::Buffer::FromVector(std::move(end_pos)), nullptr, 0, 0);

  auto df = export_pyarrow_table(
                arrow::Table::Make(std::move(schema), {chrom_data, start_data, end_data}))
                .attr("to_pandas")(nb::arg("self_destruct") = true);

  if (bin_ids.empty()) {
    return df;
  }

  auto bin_ids_data = std::make_shared<arrow::UInt64Array>(
      num_bins, arrow::Buffer::FromVector(std::move(bin_ids)), nullptr, 0, 0);

  schema =
      std::make_shared<arrow::Schema>(arrow::FieldVector{arrow::field("bin_id", arrow::uint64())});

  auto bin_ids_df =
      export_pyarrow_table(arrow::Table::Make(std::move(schema), {std::move(bin_ids_data)}))
          .attr("to_pandas")(nb::arg("self_destruct") = true);

  df.attr("index") = bin_ids_df.attr("__getitem__")("bin_id");
  return df;
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

nb::object BinTable::coord_to_bin(std::string_view chrom, std::uint32_t pos) const {
  return nb::cast(_bins->at(chrom, pos));
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
                                 const std::vector<std::uint32_t>& positions) const -> BinIDsT {
  auto np = import_module_checked("numpy");

  if (chroms.size() != positions.size()) {
    throw std::runtime_error("chroms and positions should have the same size");
  }

  auto np_array =
      np.attr("empty")(static_cast<std::int64_t>(chroms.size()), nb::arg("dtype") = "int");
  auto buffer = nb::cast<BinIDsT>(np_array);
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
  import_pyarrow_checked();
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
        return nb::make_iterator(nb::type<hictk::BinTable>(), "BinTableIterator", bins.begin(),
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

nb::object BinTable::to_df(std::string_view range, std::string_view query_type) const {
  auto pd = nb::module_::import_("pandas");

  const auto qt =
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED;
  const auto query =
      range.empty() ? hictk::GenomicInterval{}
                    : hictk::GenomicInterval::parse(_bins->chromosomes(), std::string{range}, qt);

  const auto n = compute_num_bins(*_bins, query);

  std::vector<std::uint64_t> bin_ids(n);
  std::vector<std::int32_t> chrom_ids(n);
  std::vector<std::uint32_t> starts(n);
  std::vector<std::uint32_t> ends(n);

  std::visit(
      [&](const auto& bins) {
        const auto [first_bin, last_bin] =
            range.empty() ? std::make_pair(bins.begin(), bins.end()) : bins.find_overlap(query);
        std::size_t i = 0;
        std::for_each(first_bin, last_bin, [&](const auto& bin) {
          bin_ids[i] = bin.id();
          chrom_ids[i] = static_cast<std::int32_t>(bin.chrom().id());
          starts[i] = bin.start();
          ends[i] = bin.end();
          ++i;
        });
      },
      _bins->get());

  return make_bin_table_df(chrom_names(), std::move(chrom_ids), std::move(starts), std::move(ends),
                           std::move(bin_ids));
}

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
  nb::class_<hictk::Bin>(m, "Bin", "Genomic Bin.")
      .def_prop_ro("id", [](const hictk::Bin& b) { return b.id(); })
      .def_prop_ro("rel_id", [](const hictk::Bin& b) { return b.rel_id(); })
      .def_prop_ro("chrom", [](const hictk::Bin& b) { return b.chrom().name(); })
      .def_prop_ro("start", [](const hictk::Bin& b) { return b.start(); })
      .def_prop_ro("end", [](const hictk::Bin& b) { return b.end(); })
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
         "resolution");

  bt.def("__repr__", &BinTable::repr);

  bt.def("chromosomes", &get_chromosomes_from_object<hictkpy::BinTable>,
         nb::arg("include_all") = false,
         "Get chromosomes sizes as a dictionary mapping names to sizes.");

  bt.def("bin_size", &BinTable::resolution,
         "Get the bin size for the bin table. "
         "Return 0 in case the bin table has a variable bin size.");

  bt.def("type", &BinTable::type,
         "Get the type of table underlying the BinTable object (i.e. fixed or variable).");

  bt.def("__len__", &BinTable::size, "Get the number of bins in the bin table.");

  bt.def("__iter__", &BinTable::make_iterable, nb::keep_alive<0, 1>());

  bt.def("get", &BinTable::bin_id_to_coord, nb::arg("bin_id"),
         "Get the genomic coordinate given a bin ID.");
  bt.def("get", &BinTable::bin_ids_to_coords, nb::arg("bin_ids"),
         "Get the genomic coordinates given a vector of bin IDs.");

  bt.def("get", &BinTable::coord_to_bin, nb::arg("chrom"), nb::arg("pos"),
         "Get the bin overlapping the given genomic coordinate.");
  bt.def("get", &BinTable::coords_to_bins, nb::arg("chroms"), nb::arg("pos"),
         "Get the bins overlapping the given genomic coordinates.");

  bt.def("get_id", &BinTable::coord_to_bin_id, nb::arg("chrom"), nb::arg("pos"),
         "Get the ID of the bin overlapping the given genomic coordinate.");
  bt.def("get_ids", &BinTable::coords_to_bin_ids, nb::arg("chroms"), nb::arg("pos"),
         "Get the IDs of the bins overlapping the given genomic coordinates.");

  bt.def("merge", &BinTable::merge_coords, nb::arg("df"),
         "Merge genomic coordinates corresponding to the given bin identifiers. "
         "Bin identifiers should be provided as a pandas.DataFrame with columns \"bin1_id\" and "
         "\"bin2_id\"");

  bt.def("to_df", &BinTable::to_df, nb::arg("range") = "", nb::arg("query_type") = "UCSC",
         "Convert the bin table to a pandas.DataFrame");
}

}  // namespace hictkpy
