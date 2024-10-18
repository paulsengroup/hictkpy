// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/hic_file_writer.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/reference.hpp>
#include <hictk/tmpdir.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel.hpp"
#include "hictkpy/reference.hpp"

namespace nb = nanobind;

namespace hictkpy {

HiCFileWriter::HiCFileWriter(const std::filesystem::path &path_, const ChromosomeDict &chromosomes,
                             const std::vector<std::uint32_t> &resolutions_,
                             std::string_view assembly, std::size_t n_threads,
                             std::size_t chunk_size, const std::filesystem::path &tmpdir,
                             std::uint32_t compression_lvl, bool skip_all_vs_all_matrix)
    : _tmpdir(tmpdir, true),
      _w(path_.string(), chromosome_dict_to_reference(chromosomes), resolutions_, assembly,
         n_threads, chunk_size, _tmpdir(), compression_lvl, skip_all_vs_all_matrix) {
  SPDLOG_INFO(FMT_STRING("using \"{}\" folder to store temporary file(s)"), _tmpdir());
}

HiCFileWriter::HiCFileWriter(const std::filesystem::path &path_, const ChromosomeDict &chromosomes,
                             std::uint32_t resolution, std::string_view assembly,
                             std::size_t n_threads, std::size_t chunk_size,
                             const std::filesystem::path &tmpdir, std::uint32_t compression_lvl,
                             bool skip_all_vs_all_matrix)
    : HiCFileWriter(path_, chromosomes, std::vector<std::uint32_t>{resolution}, assembly, n_threads,
                    chunk_size, tmpdir, compression_lvl, skip_all_vs_all_matrix) {}

void HiCFileWriter::finalize([[maybe_unused]] std::string_view log_lvl_str) {
  if (_finalized) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("finalize() was already called on file \"{}\""), _w.path()));
  }

#ifndef _WIN32
  // TODO fixme
  // There is something very odd going on when trying to call most spdlog functions from within
  // Python bindings on recent versions of Windows.
  // Possibly related to https://github.com/gabime/spdlog/issues/3212
  const auto log_lvl = spdlog::level::from_str(normalize_log_lvl(log_lvl_str));
  const auto previous_lvl = spdlog::default_logger()->level();
  spdlog::default_logger()->set_level(log_lvl);

  SPDLOG_INFO(FMT_STRING("finalizing file \"{}\"..."), _w.path());
#endif
  try {
    _w.serialize();
    _finalized = true;
  } catch (...) {
#ifndef _WIN32
    spdlog::default_logger()->set_level(previous_lvl);
#endif
    throw;
  }
  SPDLOG_INFO(FMT_STRING("successfully finalized \"{}\"!"), _w.path());
#ifndef _WIN32
  spdlog::default_logger()->set_level(previous_lvl);
#endif
}

std::filesystem::path HiCFileWriter::path() const noexcept {
  return std::filesystem::path{_w.path()};
}

const std::vector<std::uint32_t> &HiCFileWriter::resolutions() const noexcept {
  return _w.resolutions();
}

const hictk::Reference &HiCFileWriter::chromosomes() const { return _w.chromosomes(); }

void HiCFileWriter::add_pixels(const nb::object &df) {
  if (_finalized) {
    throw std::runtime_error(
        "caught attempt to add_pixels to a .hic file that has already been finalized!");
  }

  const auto coo_format = nb::cast<bool>(df.attr("columns").attr("__contains__")("bin1_id"));
  const auto pixels =
      coo_format ? coo_df_to_thin_pixels<float>(df, false)
                 : bg2_df_to_thin_pixels<float>(_w.bins(_w.resolutions().front()), df, false);
  SPDLOG_INFO(FMT_STRING("adding {} pixels to file \"{}\"..."), pixels.size(), _w.path());
  _w.add_pixels(_w.resolutions().front(), pixels.begin(), pixels.end());
}

std::string HiCFileWriter::repr() const {
  return fmt::format(FMT_STRING("HiCFileWriter({})"), _w.path());
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
             "Open a .hic file for writing.");

  writer.def(nb::init<const std::filesystem::path &, const ChromosomeDict &,
                      const std::vector<std::uint32_t> &, std::string_view, std::size_t,
                      std::size_t, const std::filesystem::path &, std::uint32_t, bool>(),
             nb::arg("path"), nb::arg("chromosomes"), nb::arg("resolutions"),
             nb::arg("assembly") = "unknown", nb::arg("n_threads") = 1,
             nb::arg("chunk_size") = 10'000'000,
             nb::arg("tmpdir") = hictk::internal::TmpDir::default_temp_directory_path(),
             nb::arg("compression_lvl") = 10, nb::arg("skip_all_vs_all_matrix") = false,
             "Open a .hic file for writing.");
  // NOLINTEND(*-avoid-magic-numbers)

  writer.def("__repr__", &hictkpy::HiCFileWriter::repr);

  writer.def("path", &hictkpy::HiCFileWriter::path, "Get the file path.");
  writer.def("resolutions", &hictkpy::HiCFileWriter::resolutions,
             "Get the list of resolutions in bp.");
  writer.def("chromosomes", &get_chromosomes_from_file<hictkpy::HiCFileWriter>,
             nb::arg("include_all") = false,
             "Get chromosomes sizes as a dictionary mapping names to sizes.");

  writer.def("add_pixels", &hictkpy::HiCFileWriter::add_pixels,
             nb::sig("def add_pixels(self, pixels: pd.DataFrame) -> None"), nb::arg("pixels"),
             "Add pixels from a pandas DataFrame containing pixels in COO or BG2 format (i.e. "
             "either with columns=[bin1_id, bin2_id, count] or with columns=[chrom1, start1, end1, "
             "chrom2, start2, end2, count].");
  writer.def("finalize", &hictkpy::HiCFileWriter::finalize, nb::arg("log_lvl") = "WARN",
             "Write interactions to file.");
}

}  // namespace hictkpy
