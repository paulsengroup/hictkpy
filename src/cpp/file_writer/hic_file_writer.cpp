// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/hic_file_writer.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/bin_table.hpp>
#include <hictk/file.hpp>
#include <hictk/reference.hpp>
#include <hictk/tmpdir.hpp>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/file_writer_helpers.hpp"
#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel_table.hpp"
#include "hictkpy/reference.hpp"
#include "hictkpy/table.hpp"
#include "hictkpy/to_numpy.hpp"

namespace nb = nanobind;

namespace hictkpy {

static HiCFileWriter &ctx_enter(HiCFileWriter &w) { return w; }

static void ctx_exit(HiCFileWriter &w, nb::handle exc_type, [[maybe_unused]] nb::handle exc_value,
                     [[maybe_unused]] nb::handle traceback) {
  const auto exc_raised = [exc_type]() {
    HICTKPY_GIL_SCOPED_ACQUIRE
    return !exc_type.is_none();
  }();

  if (exc_raised) {
    w.try_cleanup();
    return;
  }

  if (!w.finalized()) {
    std::ignore = w.finalize(std::nullopt);
  }
}

static ChromosomeDict get_chromosomes_checked(const hictk::BinTable &bins) {
  if (bins.type() != hictk::BinTable::Type::fixed) {
    throw std::runtime_error(
        "constructing .hic files is only supported when the BinTable has a uniform bin size");
  }

  HICTKPY_GIL_SCOPED_ACQUIRE
  ChromosomeDict chroms{};
  for (const auto &chrom : bins.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    const auto *chrom_name = chrom.name().data();
    chroms[chrom_name] = chrom.size();
  }

  return chroms;
}

static std::uint32_t get_base_resolution(const std::vector<std::uint32_t> &resolutions) {
  if (resolutions.empty()) {
    throw std::invalid_argument("please provide one or more resolutions");
  }

  return resolutions.front();
}

HiCFileWriter::HiCFileWriter(const std::filesystem::path &path_, const ChromosomeDict &chromosomes,
                             const std::vector<std::uint32_t> &resolutions_,
                             std::string_view assembly, std::size_t n_threads,
                             std::size_t chunk_size, const std::filesystem::path &tmpdir_,
                             std::uint32_t compression_lvl, bool skip_all_vs_all_matrix)
    : _path(path_.string()),
      _base_resolution(get_base_resolution(resolutions_)),
      _tmpdir(std::make_optional<hictk::internal::TmpDir>(tmpdir_, true)),
      _w(std::make_optional<hictk::hic::internal::HiCFileWriter>(
          _path, chromosome_dict_to_reference(chromosomes), resolutions_, assembly, n_threads,
          chunk_size, tmpdir(), compression_lvl, skip_all_vs_all_matrix)) {
  SPDLOG_INFO(FMT_STRING("using \"{}\" folder to store temporary file(s)"), tmpdir().string());
}

HiCFileWriter::HiCFileWriter(const std::filesystem::path &path_, const ChromosomeDict &chromosomes,
                             std::uint32_t resolution, std::string_view assembly,
                             std::size_t n_threads, std::size_t chunk_size,
                             const std::filesystem::path &tmpdir_, std::uint32_t compression_lvl,
                             bool skip_all_vs_all_matrix)
    : HiCFileWriter(path_, chromosomes, std::vector<std::uint32_t>{resolution}, assembly, n_threads,
                    chunk_size, tmpdir_, compression_lvl, skip_all_vs_all_matrix) {}

HiCFileWriter::HiCFileWriter(const std::filesystem::path &path_, const hictkpy::BinTable &bins_,
                             std::string_view assembly, std::size_t n_threads,
                             std::size_t chunk_size, const std::filesystem::path &tmpdir_,
                             std::uint32_t compression_lvl, bool skip_all_vs_all_matrix)
    : HiCFileWriter(path_, get_chromosomes_checked(*bins_.get()), bins_.get()->resolution(),
                    assembly, n_threads, chunk_size, tmpdir_, compression_lvl,
                    skip_all_vs_all_matrix) {}

File HiCFileWriter::finalize(std::optional<std::string_view> log_lvl_str) {
  if (finalized()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("finalize() was already called on file \"{}\""), _path));
  }

  if (log_lvl_str.has_value()) {
    raise_python_deprecation_warning(
        FMT_STRING("HiCFileWriter::finalize(): changing log level with argument log_lvl=\"{0}\" is "
                   "deprecated and has no effect.\n"
                   "Please use hictkpy.logging.setLevel(\"{0}\") to change the log level instead."),
        log_lvl_str);
  }
  SPDLOG_INFO(FMT_STRING("finalizing file \"{}\"..."), get().path());
  auto writer = [&]() {
    HICTKPY_GIL_SCOPED_ACQUIRE
    hictk::hic::internal::HiCFileWriter w{std::move(*_w)};
    _w.reset();
    return w;
  }();

  try {
    writer.serialize();
  } catch (...) {
    _w = std::move(writer);
    throw;
  }

  _tmpdir.reset();
  SPDLOG_INFO(FMT_STRING("successfully finalized \"{}\"!"), _path);

  return File{hictk::hic::File{_path, _base_resolution}};
}

bool HiCFileWriter::finalized() const noexcept { return !_w.has_value(); }

void HiCFileWriter::try_cleanup() noexcept {
  try {
    SPDLOG_DEBUG("HiCFileWriter::try_cleanup()");
    _w.reset();
    _tmpdir.reset();
  } catch (...) {  // NOLINT
  }
}

std::filesystem::path HiCFileWriter::path() const noexcept { return std::filesystem::path{_path}; }

auto HiCFileWriter::resolutions() const {
  return make_owning_numpy<std::int64_t>(get().resolutions());
}

const hictk::Reference &HiCFileWriter::chromosomes() const { return get().chromosomes(); }

hictkpy::BinTable HiCFileWriter::bins(std::uint32_t resolution) const {
  return hictkpy::BinTable{get().bins(resolution)};
}

void HiCFileWriter::add_pixels_from_dict(const nb::dict &columns, bool validate) {
  add_pixels(internal::make_table(columns), validate);
}

void HiCFileWriter::add_pixels_from_table(const nb::object &df, bool validate) {
  add_pixels(import_pyarrow_table(df), validate);
}

void HiCFileWriter::add_pixels(const PyArrowTable &table, bool validate) {
  if (finalized()) {
    throw std::runtime_error(
        "caught attempt to add_pixels() to a .hic file that has already been finalized!");
  }

  if (!table) {
    return;
  }

  if (table.type() != PyArrowTable::Type::BG2 && table.type() != PyArrowTable::Type::COO) {
    internal::raise_invalid_table_format();
  }

  const auto pixels = std::get<ThinPixelBuffer<float>>(
      convert_table_to_thin_pixels(get().bins(_base_resolution), table, false, float{}));

  SPDLOG_INFO(FMT_STRING("adding {} pixels to file \"{}\"..."), pixels.size(), get().path());
  get().add_pixels(_base_resolution, pixels.begin(), pixels.end(), validate);
}

std::string HiCFileWriter::repr() const {
  return fmt::format(FMT_STRING("HiCFileWriter({})"), _path);
}

const std::filesystem::path &HiCFileWriter::tmpdir() const {
  if (!_tmpdir.has_value()) {
    assert(!_w.has_value());
    throw std::runtime_error(fmt::format(
        FMT_STRING("caught an attempt to access file \"{}\", which has already been closed"),
        _path));
  }

  return (*_tmpdir)();
}

hictk::hic::internal::HiCFileWriter &HiCFileWriter::get() {
  if (!_w.has_value()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("caught an attempt to access file \"{}\", which has already been closed"),
        _path));
  }
  return *_w;
}

const hictk::hic::internal::HiCFileWriter &HiCFileWriter::get() const {
  if (!_w.has_value()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("caught an attempt to access file \"{}\", which has already been closed"),
        _path));
  }
  return *_w;
}

void HiCFileWriter::bind(nb::module_ &m) {
  auto hic = m.def_submodule("hic");

  auto writer = nb::class_<hictkpy::HiCFileWriter>(
      hic, "FileWriter", "Class representing a file handle to create .hic files.");

  // NOLINTBEGIN(*-avoid-magic-numbers)
  writer.def(nb::init<const std::filesystem::path &, const ChromosomeDict &, std::uint32_t,
                      std::string_view, std::size_t, std::size_t, const std::filesystem::path &,
                      std::uint32_t, bool>(),
             nb::arg("path"), nb::arg("chromosomes"), nb::arg("resolution"),
             nb::arg("assembly") = "unknown", nb::arg("n_threads") = 1,
             nb::arg("chunk_size") = 10'000'000,
             nb::arg("tmpdir") = hictk::internal::TmpDir::default_temp_directory_path(),
             nb::arg("compression_lvl") = 10, nb::arg("skip_all_vs_all_matrix") = false,
             "Open a .hic file for writing given a list of chromosomes with their sizes and one "
             "resolution.");

  writer.def(
      nb::init<const std::filesystem::path &, const ChromosomeDict &,
               const std::vector<std::uint32_t> &, std::string_view, std::size_t, std::size_t,
               const std::filesystem::path &, std::uint32_t, bool>(),
      nb::arg("path"), nb::arg("chromosomes"), nb::arg("resolutions"),
      nb::arg("assembly") = "unknown", nb::arg("n_threads") = 1, nb::arg("chunk_size") = 10'000'000,
      nb::arg("tmpdir") = hictk::internal::TmpDir::default_temp_directory_path(),
      nb::arg("compression_lvl") = 10, nb::arg("skip_all_vs_all_matrix") = false,
      "Open a .hic file for writing given a list of chromosomes with their sizes and one or more "
      "resolutions.");

  writer.def(
      nb::init<const std::filesystem::path &, const hictkpy::BinTable &, std::string_view,
               std::size_t, std::size_t, const std::filesystem::path &, std::uint32_t, bool>(),
      nb::arg("path"), nb::arg("bins"), nb::arg("assembly") = "unknown", nb::arg("n_threads") = 1,
      nb::arg("chunk_size") = 10'000'000,
      nb::arg("tmpdir") = hictk::internal::TmpDir::default_temp_directory_path(),
      nb::arg("compression_lvl") = 10, nb::arg("skip_all_vs_all_matrix") = false,
      "Open a .hic file for writing given a BinTable. Only BinTable with a fixed bin size are "
      "supported.");
  // NOLINTEND(*-avoid-magic-numbers)

  writer.def("__repr__", &hictkpy::HiCFileWriter::repr, nb::rv_policy::move);

  writer.def("__enter__", &ctx_enter, nb::rv_policy::reference_internal);

  writer.def("__exit__", &ctx_exit,
             // clang-format off
             nb::call_guard<nb::gil_scoped_release>(),
             nb::arg("exc_type") = nb::none(),
             nb::arg("exc_value") = nb::none(),
             nb::arg("traceback") = nb::none()
             // clang-format on
  );

  writer.def("path", &hictkpy::HiCFileWriter::path, "Get the file path.", nb::rv_policy::move);
  writer.def("resolutions", &hictkpy::HiCFileWriter::resolutions,
             "Get the list of resolutions in bp.", nb::rv_policy::take_ownership);
  writer.def("chromosomes", &get_chromosomes_from_object<hictkpy::HiCFileWriter>,
             nb::arg("include_ALL") = false,
             "Get the chromosome sizes as a dictionary mapping names to sizes.",
             nb::rv_policy::take_ownership);
  writer.def("bins", &hictkpy::HiCFileWriter::bins, "Get table of bins for the given resolution.",
             nb::sig("def bins(self, resolution: int) -> hictkpy.BinTable"), nb::rv_policy::move);

  writer.def(
      "add_pixels", &hictkpy::HiCFileWriter::add_pixels_from_table,
      nb::call_guard<nb::gil_scoped_release>(),
      nb::sig("def add_pixels(self, pixels: pandas.DataFrame | pyarrow.Table, validate: bool = "
              "True) -> None"),
      nb::arg("pixels"), nb::arg("validate") = true,
      "Add pixels from a pandas.DataFrame or pyarrow.Table containing pixels in COO or BG2 format "
      "(i.e. either with columns=[bin1_id, bin2_id, count] or with "
      "columns=[chrom1, start1, end1, chrom2, start2, end2, count]).\n"
      "When sorted is True, pixels are assumed to be sorted by their genomic coordinates in "
      "ascending order.\n"
      "When validate is True, hictkpy will perform some basic sanity checks on the given "
      "pixels before adding them to the .hic file.");
  writer.def("add_pixels_from_dict", &hictkpy::HiCFileWriter::add_pixels_from_dict,
             nb::call_guard<nb::gil_scoped_release>(),
             nb::sig("def add_pixels_from_dict(self, columns: Dict[str, Iterable[str | int | "
                     "float]], validate: bool = True) -> None"),
             nb::arg("pixels"), nb::arg("validate") = true,
             "Add pixels from a dictionary containing containing columns corresponding to pixels "
             "in COO or BG2 format (i.e. either with keys=[bin1_id, bin2_id, count] or with "
             "keys=[chrom1, start1, end1, chrom2, start2, end2, count]).\n"
             "When sorted is True, pixels are assumed to be sorted by their genomic coordinates in "
             "ascending order.\n"
             "When validate is True, hictkpy will perform some basic sanity checks on the given "
             "pixels before adding them to the Cooler file.");
  writer.def("finalize", &hictkpy::HiCFileWriter::finalize,
             nb::call_guard<nb::gil_scoped_release>(), nb::arg("log_lvl") = nb::none(),
             "Write interactions to file.", nb::rv_policy::move);
}

}  // namespace hictkpy
