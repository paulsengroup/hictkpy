// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/cooler_file_writer.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <hictk/cooler/cooler.hpp>
#include <hictk/reference.hpp>
#include <hictk/type_traits.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel.hpp"
#include "hictkpy/reference.hpp"

namespace nb = nanobind;

namespace hictkpy {

CoolerFileWriter::CoolerFileWriter(std::string_view path_, const nb::dict &chromosomes_,
                                   std::uint32_t resolution_, std::string_view assembly,
                                   const std::filesystem::path &tmpdir,
                                   std::uint32_t compression_lvl)
    : _path(std::string{path_}),
      _tmpdir(tmpdir / (_path + ".tmp")),
      _w(create_file(_path, chromosome_dict_to_reference(chromosomes_), resolution_, assembly,
                     _tmpdir())),
      _compression_lvl(compression_lvl) {
  if (std::filesystem::exists(_path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to create .cool file \"{}\": file already exists"), path()));
  }
}

std::string_view CoolerFileWriter::path() const noexcept { return _path; }

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

void CoolerFileWriter::add_pixels(const nb::object &df) {
  if (!_w.has_value()) {
    throw std::runtime_error(
        "caught attempt to add_pixels to a .cool file that has already been finalized!");
  }

  const auto coo_format = nb::cast<bool>(df.attr("columns").attr("__contains__")("bin1_id"));

  const auto cell_id = fmt::to_string(_w->cells().size());
  auto attrs = hictk::cooler::Attributes::init(_w->resolution());
  attrs.assembly = _w->attributes().assembly;

  const auto dtype = df.attr("__getitem__")("count").attr("dtype");
  const auto dtype_str = nb::cast<std::string>(dtype.attr("__str__")());
  const auto var = map_dtype_to_type(dtype_str);

  std::visit(
      [&](const auto &n) {
        using N = hictk::remove_cvref_t<decltype(n)>;
        const auto pixels = coo_format ? coo_df_to_thin_pixels<N>(df, true)
                                       : bg2_df_to_thin_pixels<N>(_w->bins(), df, true);

        auto clr = _w->create_cell<N>(cell_id, std::move(attrs),
                                      hictk::cooler::DEFAULT_HDF5_CACHE_SIZE * 4, 1);
        clr.append_pixels(pixels.begin(), pixels.end());

        clr.flush();
      },
      var);
}

void CoolerFileWriter::finalize([[maybe_unused]] std::string_view log_lvl_str,
                                std::size_t chunk_size, std::size_t update_freq) {
  if (_finalized) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("finalize() was already called on file \"{}\""), _path));
  }

  if (chunk_size == 0) {
    throw std::runtime_error("chunk_size must be greater than 0");
  }

  assert(_w.has_value());
  // NOLINTBEGIN(*-unchecked-optional-access)
#ifndef _WIN32
  // TODO fixme
  // There is something very odd going on when trying to call most spdlog functions from within
  // Python bindings on recent versions of Windows.
  // Possibly related to https://github.com/gabime/spdlog/issues/3212
  const auto log_lvl = spdlog::level::from_str(normalize_log_lvl(log_lvl_str));
  const auto previous_lvl = spdlog::default_logger()->level();
  spdlog::default_logger()->set_level(log_lvl);
#endif
  try {
    std::visit(
        [&](const auto &num) {
          using N = hictk::remove_cvref_t<decltype(num)>;
          _w->aggregate<N>(_path, false, _compression_lvl, chunk_size, update_freq);
        },
        _w->open("0").pixel_variant());
  } catch (...) {
#ifndef _WIN32
    spdlog::default_logger()->set_level(previous_lvl);
#endif
    throw;
  }

  _finalized = true;
#ifndef _WIN32
  spdlog::default_logger()->set_level(previous_lvl);
#endif
  const std::string sclr_path{_w->path()};
  _w.reset();
  std::filesystem::remove(sclr_path);  // NOLINT
  // NOLINTEND(*-unchecked-optional-access)
}

hictk::cooler::SingleCellFile CoolerFileWriter::create_file(std::string_view path,
                                                            const hictk::Reference &chromosomes,
                                                            std::uint32_t resolution,
                                                            std::string_view assembly,
                                                            const std::filesystem::path &tmpdir) {
  auto attrs = hictk::cooler::SingleCellAttributes::init(resolution);
  attrs.assembly = assembly;
  return hictk::cooler::SingleCellFile::create(tmpdir / std::filesystem::path{path}.filename(),
                                               chromosomes, resolution, false, std::move(attrs));
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
  writer.def(nb::init<std::string_view, nb::dict, std::uint32_t, std::string_view,
                      const std::filesystem::path &, std::uint32_t>(),
             nb::arg("path"), nb::arg("chromosomes"), nb::arg("resolution"),
             nb::arg("assembly") = "unknown",
             nb::arg("tmpdir") = std::filesystem::temp_directory_path().string(),
             nb::arg("compression_lvl") = 6, "Open a .cool file for writing.");
  // NOLINTEND(*-avoid-magic-numbers)

  writer.def("__repr__", &hictkpy::CoolerFileWriter::repr);

  writer.def("path", &hictkpy::CoolerFileWriter::path, "Get the file path.");
  writer.def("resolutions", &hictkpy::CoolerFileWriter::resolution, "Get the resolution in bp.");
  writer.def("chromosomes", &get_chromosomes_from_file<hictkpy::CoolerFileWriter>,
             nb::arg("include_all") = false,
             "Get chromosomes sizes as a dictionary mapping names to sizes.");

  writer.def("add_pixels", &hictkpy::CoolerFileWriter::add_pixels,
             nb::sig("def add_pixels(self, pixels: pandas.DataFrame)"), nb::arg("pixels"),
             "Add pixels from a pandas DataFrame containing pixels in COO or BG2 format (i.e. "
             "either with columns=[bin1_id, bin2_id, count] or with columns=[chrom1, start1, end1, "
             "chrom2, start2, end2, count].");
  // NOLINTBEGIN(*-avoid-magic-numbers)
  writer.def("finalize", &hictkpy::CoolerFileWriter::finalize, nb::arg("log_lvl") = "WARN",
             nb::arg("chunk_size") = 500'000, nb::arg("update_frequency") = 10'000'000,
             "Write interactions to file.");
  // NOLINTEND(*-avoid-magic-numbers)
}
}  // namespace hictkpy
