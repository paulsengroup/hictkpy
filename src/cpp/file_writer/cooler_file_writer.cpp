// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/cooler_file_writer.hpp"

#include <arrow/table.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/bin_table.hpp>
#include <hictk/cooler/cooler.hpp>
#include <hictk/file.hpp>
#include <hictk/reference.hpp>
#include <hictk/tmpdir.hpp>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <variant>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/file.hpp"
#include "hictkpy/file_writer_helpers.hpp"
#include "hictkpy/locking.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel_table.hpp"
#include "hictkpy/reference.hpp"
#include "hictkpy/table.hpp"
#include "hictkpy/type.hpp"
#include "hictkpy/variant.hpp"

namespace nb = nanobind;

namespace hictkpy {

static CoolerFileWriter &ctx_enter(CoolerFileWriter &w) { return w; }

static void ctx_exit(CoolerFileWriter &w, nb::handle exc_type,
                     [[maybe_unused]] nb::handle exc_value, [[maybe_unused]] nb::handle traceback) {
  const auto exc_raised = [exc_type]() {
    HICTKPY_GIL_SCOPED_ACQUIRE
    return !exc_type.is_none();
  }();

  if (exc_raised) {
    w.try_cleanup();
    return;
  }

  if (!w.finalized()) {  // NOLINTNEXTLINE(*-avoid-magic-numbers)
    std::ignore = w.finalize(std::nullopt, 500'000, 10'000'000);
  }
}

CoolerFileWriter::CoolerFileWriter(std::filesystem::path path_, const hictkpy::BinTable &bins_,
                                   std::string_view assembly, const std::filesystem::path &tmpdir_,
                                   std::uint32_t compression_lvl)
    : _path(std::move(path_)),
      _tmpdir(std::make_optional<hictk::internal::TmpDir>(tmpdir_, true)),
      _w(create_file(_path.string(), *bins_.get(), assembly, tmpdir())),
      _compression_lvl(compression_lvl) {
  if (std::filesystem::exists(_path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to create .cool file \"{}\": file already exists"), path()));
  }

  SPDLOG_INFO(FMT_STRING("using \"{}\" folder to store temporary file(s)"), tmpdir().string());
}

CoolerFileWriter::CoolerFileWriter(std::filesystem::path path_, const ChromosomeDict &chromosomes_,
                                   std::uint32_t resolution_, std::string_view assembly,
                                   const std::filesystem::path &tmpdir,
                                   std::uint32_t compression_lvl)
    : CoolerFileWriter(std::move(path_), BinTable{chromosomes_, resolution_}, assembly, tmpdir,
                       compression_lvl) {}

const std::filesystem::path &CoolerFileWriter::path() const noexcept { return _path; }

std::uint32_t CoolerFileWriter::resolution() const { return get().resolution(); }

const hictk::Reference &CoolerFileWriter::chromosomes() const { return get().chromosomes(); }

std::shared_ptr<const hictk::BinTable> CoolerFileWriter::bins_ptr() const {
  return get().bins_ptr();
}

void CoolerFileWriter::add_pixels_from_dict(const nb::dict &columns, bool sorted, bool validate) {
  add_pixels(internal::make_table(columns), sorted, validate);
}

void CoolerFileWriter::add_pixels_from_table(const nb::object &df, bool sorted, bool validate) {
  add_pixels(import_pyarrow_table(df), sorted, validate);
}

void CoolerFileWriter::add_pixels(const PyArrowTable &table, bool sorted, bool validate) {
  if (finalized()) {
    throw std::runtime_error(
        "caught attempt to add_pixels() to a .cool file that has already been finalized!");
  }

  if (!table) {
    return;
  }

  if (table.type() != PyArrowTable::Type::BG2 && table.type() != PyArrowTable::Type::COO) {
    internal::raise_invalid_table_format();
  }

  const auto count_type = internal::infer_count_type(table.get());
  const auto pixel_buff = convert_table_to_thin_pixels(get().bins(), table, !sorted, count_type);

  const auto cell_id = fmt::to_string(get().cells().size());
  auto attrs = hictk::cooler::Attributes::init(resolution());
  attrs.assembly = get().attributes().assembly;

  std::visit(
      [&](const auto &pixels) {
        using N = remove_cvref_t<decltype(pixels.front().count)>;
        HICTKPY_LOCK_COOLER_MTX_SCOPED
        auto clr = [&]() {
          return get().create_cell<N>(cell_id, std::move(attrs),
                                      hictk::cooler::DEFAULT_HDF5_CACHE_SIZE * 4, 1);
        }();

        SPDLOG_INFO(FMT_STRING("adding {} pixels of type {} to file \"{}\"..."), pixels.size(),
                    type_to_str<N>(), clr.uri());

        clr.append_pixels(pixels.begin(), pixels.end(), validate);
        clr.flush();
      },
      pixel_buff);
}

File CoolerFileWriter::finalize(std::optional<std::string_view> log_lvl_str, std::size_t chunk_size,
                                std::size_t update_freq) {
  if (finalized()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("finalize() was already called on file \"{}\""), _path.string()));
  }

  if (chunk_size == 0) {
    throw std::runtime_error("chunk_size must be greater than 0");
  }

  if (log_lvl_str.has_value()) {
    raise_python_deprecation_warning(
        FMT_STRING(
            "CoolerFileWriter::finalize(): changing log level with argument log_lvl=\"{0}\" is "
            "deprecated and has no effect.\n"
            "Please use hictkpy.logging.setLevel(\"{0}\") to change the log level instead."),
        log_lvl_str);
  }

  SPDLOG_INFO(FMT_STRING("finalizing file \"{}\"..."), _path.string());
  NumericDtype count_type{std::int32_t{}};
  auto writer = [&]() {
    HICTKPY_GIL_SCOPED_ACQUIRE
    HICTKPY_LOCK_COOLER_MTX_SCOPED
    if (!get().cells().empty()) {
      std::visit(
          [&](const auto &t) {
            using N = remove_cvref_t<decltype(t)>;
            if constexpr (std::is_same_v<N, long double>) {
              count_type = double{};
            } else {
              count_type = N{};
            }
          },
          get().open("0").pixel_variant());
    }

    decltype(_w) w = std::move(_w);
    _w.reset();
    return w;
  }();

  auto clr = std::visit(
      [&]([[maybe_unused]] const auto &num) {
        try {
          using N = remove_cvref_t<decltype(num)>;

          HICTKPY_LOCK_COOLER_MTX_SCOPED
          SPDLOG_DEBUG(FMT_STRING("aggregating file \"{}\" and writing results to file \"{}\"..."),
                       writer->path(), _path.string());

          // NOLINTNEXTLINE(*-unchecked-optional-access)
          return writer->aggregate<N>(_path.string(), false, _compression_lvl, chunk_size,
                                      update_freq);
        } catch (...) {
          _w = std::move(writer);
          [[maybe_unused]] std::error_code ec{};
          std::filesystem::remove(_path, ec);  // NOLINT
          throw;
        }
      },
      count_type);

  // NOLINTNEXTLINE(*-unchecked-optional-access)
  SPDLOG_INFO(FMT_STRING("merged {} cooler(s) into file \"{}\""), writer->cells().size(), _path);

  {
    HICTKPY_GIL_SCOPED_ACQUIRE
    _w = std::move(writer);
    reset();
  }

  HICTKPY_LOCK_COOLER_MTX_SCOPED
  return File{std::move(clr)};
}

bool CoolerFileWriter::finalized() const noexcept { return !_w.has_value(); }

void CoolerFileWriter::try_cleanup() noexcept {
  try {
    SPDLOG_DEBUG("CoolerFileWriter::try_cleanup()");
    reset();
  } catch (...) {  // NOLINT
  }
}

std::optional<hictk::cooler::SingleCellFile> CoolerFileWriter::create_file(
    std::string_view path, const hictk::BinTable &bins, std::string_view assembly,
    const std::filesystem::path &tmpdir) {
  using namespace hictk::cooler;
  auto attrs = SingleCellAttributes::init(bins.resolution());
  attrs.assembly = assembly;

  HICTKPY_LOCK_COOLER_MTX_SCOPED
  return std::make_optional<SingleCellFile>(SingleCellFile::create(
      tmpdir / std::filesystem::path{path}.filename(), bins, false, std::move(attrs)));
}

void CoolerFileWriter::reset() {
  const auto tmpfile = [&]() {
    HICTKPY_LOCK_COOLER_MTX_SCOPED
    auto tmpfile_ = get().path();
    _w.reset();
    return tmpfile_;
  }();

  _tmpdir.reset();
  std::filesystem::remove(tmpfile);  // NOLINT
}

const std::filesystem::path &CoolerFileWriter::tmpdir() const {
  if (!_tmpdir.has_value()) {
    assert(!_w.has_value());
    throw std::runtime_error(fmt::format(
        FMT_STRING("caught an attempt to access file \"{}\", which has already been closed"),
        _path));
  }

  return (*_tmpdir)();
}

hictk::cooler::SingleCellFile &CoolerFileWriter::get() {
  if (!_w.has_value()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("caught an attempt to access file \"{}\", which has already been closed"),
        _path));
  }
  return *_w;
}

const hictk::cooler::SingleCellFile &CoolerFileWriter::get() const {
  if (!_w.has_value()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("caught an attempt to access file \"{}\", which has already been closed"),
        _path));
  }
  return *_w;
}

std::string CoolerFileWriter::repr() const {
  return fmt::format(FMT_STRING("CoolFileWriter({})"), _path.string());
}

void CoolerFileWriter::bind(nb::module_ &m) {
  auto cooler = m.def_submodule("cooler");

  auto writer = nb::class_<hictkpy::CoolerFileWriter>(
      cooler, "FileWriter", "Class representing a file handle to create .cool files.");

  // NOLINTBEGIN(*-avoid-magic-numbers)
  writer.def(nb::init<std::filesystem::path, const ChromosomeDict &, std::uint32_t,
                      std::string_view, const std::filesystem::path &, std::uint32_t>(),
             nb::arg("path"), nb::arg("chromosomes"), nb::arg("resolution"),
             nb::arg("assembly") = "unknown",
             nb::arg("tmpdir") = hictk::internal::TmpDir::default_temp_directory_path(),
             nb::arg("compression_lvl") = 6,
             "Open a .cool file for writing given a list of chromosomes with their sizes and a "
             "resolution.");
  writer.def(nb::init<std::filesystem::path, const hictkpy::BinTable &, std::string_view,
                      const std::filesystem::path &, std::uint32_t>(),
             nb::arg("path"), nb::arg("bins"), nb::arg("assembly") = "unknown",
             nb::arg("tmpdir") = hictk::internal::TmpDir::default_temp_directory_path(),
             nb::arg("compression_lvl") = 6,
             "Open a .cool file for writing given a table of bins.");
  // NOLINTEND(*-avoid-magic-numbers)

  writer.def("__repr__", &hictkpy::CoolerFileWriter::repr, nb::rv_policy::move);

  writer.def("__enter__", &ctx_enter, nb::rv_policy::reference_internal);

  writer.def("__exit__", &ctx_exit,
             // clang-format off
             nb::call_guard<nb::gil_scoped_release>(),
             nb::arg("exc_type") = nb::none(),
             nb::arg("exc_value") = nb::none(),
             nb::arg("traceback") = nb::none()
             // clang-format on
  );

  writer.def("path", &hictkpy::CoolerFileWriter::path, "Get the file path.", nb::rv_policy::copy);
  writer.def("resolution", &hictkpy::CoolerFileWriter::resolution, "Get the resolution in bp.");
  writer.def("chromosomes", &get_chromosomes_from_object<hictkpy::CoolerFileWriter>,
             nb::arg("include_ALL") = false,
             "Get the chromosome sizes as a dictionary mapping names to sizes.",
             nb::rv_policy::take_ownership);
  writer.def("bins", &get_bins_from_object<hictkpy::CoolerFileWriter>, "Get table of bins.",
             nb::sig("def bins(self) -> hictkpy.BinTable"), nb::rv_policy::move);

  writer.def(
      "add_pixels", &hictkpy::CoolerFileWriter::add_pixels_from_table,
      nb::call_guard<nb::gil_scoped_release>(),
      nb::sig(
          "def add_pixels(self, pixels: pandas.DataFrame | pyarrow.Table, sorted: bool = False, "
          "validate: bool = True) -> None"),
      nb::arg("pixels"), nb::arg("sorted") = false, nb::arg("validate") = true,
      "Add pixels from a pandas.DataFrame or pyarrow.Table containing pixels in COO or BG2 format "
      "(i.e. either with columns=[bin1_id, bin2_id, count] or with "
      "columns=[chrom1, start1, end1, chrom2, start2, end2, count]).\n"
      "When sorted is True, pixels are assumed to be sorted by their genomic coordinates in "
      "ascending order.\n"
      "When validate is True, hictkpy will perform some basic sanity checks on the given "
      "pixels before adding them to the Cooler file.");
  writer.def("add_pixels_from_dict", &hictkpy::CoolerFileWriter::add_pixels_from_dict,
             nb::call_guard<nb::gil_scoped_release>(),
             nb::sig("def add_pixels_from_dict(self, columns: Dict[str, Iterable[str | int | "
                     "float]], sorted: bool = False, "
                     "validate: bool = True) -> None"),
             nb::arg("pixels"), nb::arg("sorted") = false, nb::arg("validate") = true,
             "Add pixels from a dictionary containing containing columns corresponding to pixels "
             "in COO or BG2 format (i.e. either with keys=[bin1_id, bin2_id, count] or with "
             "keys=[chrom1, start1, end1, chrom2, start2, end2, count]).\n"
             "When sorted is True, pixels are assumed to be sorted by their genomic coordinates in "
             "ascending order.\n"
             "When validate is True, hictkpy will perform some basic sanity checks on the given "
             "pixels before adding them to the Cooler file.");
  // NOLINTBEGIN(*-avoid-magic-numbers)
  writer.def("finalize", &hictkpy::CoolerFileWriter::finalize,
             nb::call_guard<nb::gil_scoped_release>(), nb::arg("log_lvl") = nb::none(),
             nb::arg("chunk_size") = 500'000, nb::arg("update_frequency") = 10'000'000,
             "Write interactions to file.", nb::rv_policy::move);
  // NOLINTEND(*-avoid-magic-numbers)
}
}  // namespace hictkpy
