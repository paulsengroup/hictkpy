// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/cooler_file_writer.hpp"

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <hictk/cooler/cooler.hpp>
#include <hictk/file.hpp>
#include <hictk/reference.hpp>
#include <hictk/tmpdir.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/cooler_mtx.hpp"
#include "hictkpy/file.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel.hpp"
#include "hictkpy/reference.hpp"
#include "hictkpy/type.hpp"

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

  if (!w.finalized()) {
    std::ignore = w.finalize("WARN", 500'000, 10'000'000);
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

void CoolerFileWriter::add_pixels(const nb::object &df, bool sorted, bool validate) {
  if (finalized()) {
    throw std::runtime_error(
        "caught attempt to add_pixels() to a .cool file that has already been finalized!");
  }

  const auto cell_id = fmt::to_string(get().cells().size());
  auto attrs = hictk::cooler::Attributes::init(resolution());
  attrs.assembly = get().attributes().assembly;

  const auto coo_format = [&]() {
    HICTKPY_GIL_SCOPED_ACQUIRE
    return nb::cast<bool>(df.attr("columns").attr("__contains__")("bin1_id"));
  }();

  const auto dtype_str = [&]() {
    HICTKPY_GIL_SCOPED_ACQUIRE
    auto dtype = df.attr("__getitem__")("count").attr("dtype");
    return nb::cast<std::string>(dtype.attr("__str__")());
  }();

  const auto var = map_py_numeric_to_cpp_type(dtype_str);

  std::visit(
      [&](const auto &n) {
        using N = remove_cvref_t<decltype(n)>;
        const auto pixels = coo_format ? coo_df_to_thin_pixels<N>(df, !sorted)
                                       : bg2_df_to_thin_pixels<N>(_w->bins(), df, !sorted);

        auto clr = [&]() {
          HICTKPY_LOCK_COOLER_MTX_SCOPED
          return get().create_cell<N>(cell_id, std::move(attrs),
                                      hictk::cooler::DEFAULT_HDF5_CACHE_SIZE * 4, 1);
        }();

        SPDLOG_INFO(FMT_STRING("adding {} pixels of type {} to file \"{}\"..."), pixels.size(),
                    dtype_str, clr.uri());
        {
          HICTKPY_LOCK_COOLER_MTX_SCOPED
          clr.append_pixels(pixels.begin(), pixels.end(), validate);
          clr.flush();
        }
      },
      var);
}

File CoolerFileWriter::finalize(std::string_view log_lvl_str, std::size_t chunk_size,
                                std::size_t update_freq) {
  if (finalized()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("finalize() was already called on file \"{}\""), _path.string()));
  }

  if (chunk_size == 0) {
    throw std::runtime_error("chunk_size must be greater than 0");
  }

  // TODO changing log levels in this way is problematic. Need to keep 1 logger per thread?
  // const auto log_lvl = spdlog::level::from_str(normalize_log_lvl(log_lvl_str));
  // const auto previous_lvl = spdlog::default_logger()->level();
  // spdlog::default_logger()->set_level(log_lvl);

  SPDLOG_INFO(FMT_STRING("finalizing file \"{}\"..."), _path.string());
  hictk::internal::NumericVariant count_type{std::int32_t{}};
  if (!get().cells().empty()) {
    HICTKPY_LOCK_COOLER_MTX_SCOPED
    count_type = get().open("0").pixel_variant();
  }

  try {
    std::visit(
        [&](const auto &num) {
          using N = remove_cvref_t<decltype(num)>;
          SPDLOG_DEBUG(FMT_STRING("aggregating file \"{}\" and writing results to file \"{}\"..."),
                       get().path(), _path.string());
          HICTKPY_LOCK_COOLER_MTX_SCOPED
          std::ignore =  // NOLINTNEXTLINE(*-unchecked-optional-access)
              get().aggregate<N>(_path.string(), false, _compression_lvl, chunk_size, update_freq);
        },
        count_type);

  } catch (...) {
    // spdlog::default_logger()->set_level(previous_lvl);
    throw;
  }

  SPDLOG_INFO(FMT_STRING("merged {} cooler(s) into file \"{}\""), get().cells().size(), _path);
  // spdlog::default_logger()->set_level(previous_lvl);

  reset();

  HICTKPY_LOCK_COOLER_MTX_SCOPED
  hictk::cooler::File clr{_path.string()};

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
  const std::string tmpfile{get().path()};
  {
    HICTKPY_LOCK_COOLER_MTX_SCOPED
    _w.reset();
  }
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
             nb::call_guard<nb::gil_scoped_release>(), nb::arg("path"), nb::arg("chromosomes"),
             nb::arg("resolution"), nb::arg("assembly") = "unknown",
             nb::arg("tmpdir") = hictk::internal::TmpDir::default_temp_directory_path(),
             nb::arg("compression_lvl") = 6,
             "Open a .cool file for writing given a list of chromosomes with their sizes and a "
             "resolution.");
  writer.def(nb::init<std::filesystem::path, const hictkpy::BinTable &, std::string_view,
                      const std::filesystem::path &, std::uint32_t>(),
             nb::call_guard<nb::gil_scoped_release>(), nb::arg("path"), nb::arg("bins"),
             nb::arg("assembly") = "unknown",
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

  writer.def("add_pixels", &hictkpy::CoolerFileWriter::add_pixels,
             nb::call_guard<nb::gil_scoped_release>(),
             nb::sig("def add_pixels(self, pixels: pandas.DataFrame, sorted: bool = False, "
                     "validate: bool = True) -> None"),
             nb::arg("pixels"), nb::arg("sorted") = false, nb::arg("validate") = true,
             "Add pixels from a pandas DataFrame containing pixels in COO or BG2 format (i.e. "
             "either with columns=[bin1_id, bin2_id, count] or with columns=[chrom1, start1, end1, "
             "chrom2, start2, end2, count].\n"
             "When sorted is True, pixels are assumed to be sorted by their genomic coordinates in "
             "ascending order.\n"
             "When validate is True, hictkpy will perform some basic sanity checks on the given "
             "pixels before adding them to the Cooler file.");
  // NOLINTBEGIN(*-avoid-magic-numbers)
  writer.def("finalize", &hictkpy::CoolerFileWriter::finalize,
             nb::call_guard<nb::gil_scoped_release>(), nb::arg("log_lvl") = "WARN",
             nb::arg("chunk_size") = 500'000, nb::arg("update_frequency") = 10'000'000,
             "Write interactions to file.", nb::rv_policy::move);
  // NOLINTEND(*-avoid-magic-numbers)
}
}  // namespace hictkpy
