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
#include <hictk/type_traits.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include "hictkpy/bin_table.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel.hpp"
#include "hictkpy/reference.hpp"

namespace nb = nanobind;

namespace hictkpy {

CoolerFileWriter::CoolerFileWriter(std::filesystem::path path_, const hictkpy::BinTable &bins_,
                                   std::string_view assembly, const std::filesystem::path &tmpdir,
                                   std::uint32_t compression_lvl)
    : _path(std::move(path_)),
      _tmpdir(tmpdir, true),
      _w(create_file(_path.string(), *bins_.get(), assembly, _tmpdir())),
      _compression_lvl(compression_lvl) {
  if (std::filesystem::exists(_path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to create .cool file \"{}\": file already exists"), path()));
  }

  SPDLOG_INFO(FMT_STRING("using \"{}\" folder to store temporary file(s)"), _tmpdir());
}

CoolerFileWriter::CoolerFileWriter(std::filesystem::path path_, const ChromosomeDict &chromosomes_,
                                   std::uint32_t resolution_, std::string_view assembly,
                                   const std::filesystem::path &tmpdir,
                                   std::uint32_t compression_lvl)
    : CoolerFileWriter(std::move(path_), BinTable{chromosomes_, resolution_}, assembly, tmpdir,
                       compression_lvl) {}

const std::filesystem::path &CoolerFileWriter::path() const noexcept { return _path; }

std::uint32_t CoolerFileWriter::resolution() const noexcept {
  if (_w.has_value()) {
    return _w->resolution();
  }
  return 0;
}

const hictk::Reference &CoolerFileWriter::chromosomes() const {
  if (_w.has_value()) {
    return _w->chromosomes();
  }

  const static hictk::Reference ref{};
  return ref;
}

std::shared_ptr<const hictk::BinTable> CoolerFileWriter::bins_ptr() const noexcept {
  if (!_w) {
    return {};
  }

  return _w->bins_ptr();
}

void CoolerFileWriter::add_pixels(const nb::object &df, bool sorted, bool validate) {
  if (!_w.has_value()) {
    throw std::runtime_error(
        "caught attempt to add_pixels to a .cool file that has already been finalized!");
  }

  const auto cell_id = fmt::to_string(_w->cells().size());
  auto attrs = hictk::cooler::Attributes::init(_w->resolution());
  attrs.assembly = _w->attributes().assembly;

  auto lck = std::make_optional<nb::gil_scoped_acquire>();
  const auto coo_format = nb::cast<bool>(df.attr("columns").attr("__contains__")("bin1_id"));

  const auto dtype = df.attr("__getitem__")("count").attr("dtype");
  const auto dtype_str = nb::cast<std::string>(dtype.attr("__str__")());
  const auto var = map_dtype_to_type(dtype_str);

  std::visit(
      [&](const auto &n) {
        using N = remove_cvref_t<decltype(n)>;
        const auto pixels = coo_format ? coo_df_to_thin_pixels<N>(df, !sorted)
                                       : bg2_df_to_thin_pixels<N>(_w->bins(), df, !sorted);
        lck.reset();

        auto clr = _w->create_cell<N>(cell_id, std::move(attrs),
                                      hictk::cooler::DEFAULT_HDF5_CACHE_SIZE * 4, 1);

        SPDLOG_INFO(FMT_STRING("adding {} pixels of type {} to file \"{}\"..."), pixels.size(),
                    dtype_str, clr.uri());
        clr.append_pixels(pixels.begin(), pixels.end(), validate);

        clr.flush();
      },
      var);
}

hictk::File CoolerFileWriter::finalize(std::string_view log_lvl_str, std::size_t chunk_size,
                                       std::size_t update_freq) {
  if (_finalized) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("finalize() was already called on file \"{}\""), _path));
  }

  if (chunk_size == 0) {
    throw std::runtime_error("chunk_size must be greater than 0");
  }

  assert(_w.has_value());
  // NOLINTBEGIN(*-unchecked-optional-access)
  const auto log_lvl = spdlog::level::from_str(normalize_log_lvl(log_lvl_str));
  const auto previous_lvl = spdlog::default_logger()->level();
  spdlog::default_logger()->set_level(log_lvl);

  SPDLOG_INFO(FMT_STRING("finalizing file \"{}\"..."), _path);
  hictk::internal::NumericVariant count_type{};
  if (_w->cells().empty()) {
    count_type = std::int32_t{};
  } else {
    count_type = _w->open("0").pixel_variant();
  }
  try {
    std::visit(
        [&](const auto &num) {
          using N = remove_cvref_t<decltype(num)>;
          _w->aggregate<N>(_path.string(), false, _compression_lvl, chunk_size, update_freq);
        },
        count_type);
  } catch (...) {
    spdlog::default_logger()->set_level(previous_lvl);
    throw;
  }

  _finalized = true;
  SPDLOG_INFO(FMT_STRING("merged {} cooler(s) into file \"{}\""), _w->cells().size(), _path);
  spdlog::default_logger()->set_level(previous_lvl);

  const std::string sclr_path{_w->path()};
  _w.reset();
  std::filesystem::remove(sclr_path);  // NOLINT
  // NOLINTEND(*-unchecked-optional-access)

  return hictk::File{_path.string()};
}

hictk::cooler::SingleCellFile CoolerFileWriter::create_file(std::string_view path,
                                                            const hictk::BinTable &bins,
                                                            std::string_view assembly,
                                                            const std::filesystem::path &tmpdir) {
  auto attrs = hictk::cooler::SingleCellAttributes::init(bins.resolution());
  attrs.assembly = assembly;
  return hictk::cooler::SingleCellFile::create(tmpdir / std::filesystem::path{path}.filename(),
                                               bins, false, std::move(attrs));
}

std::string CoolerFileWriter::repr() const {
  if (!_w.has_value()) {
    return "CoolFileWriter()";
  }
  return fmt::format(FMT_STRING("CoolFileWriter({})"), _w->path());
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

  writer.def("path", &hictkpy::CoolerFileWriter::path, "Get the file path.", nb::rv_policy::copy);
  writer.def("resolution", &hictkpy::CoolerFileWriter::resolution, "Get the resolution in bp.");
  writer.def("chromosomes", &get_chromosomes_from_object<hictkpy::CoolerFileWriter>,
             nb::arg("include_ALL") = false,
             "Get chromosomes sizes as a dictionary mapping names to sizes.",
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
