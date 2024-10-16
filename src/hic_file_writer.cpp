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
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictkpy/nanobind.hpp"
#include "hictkpy/pixel.hpp"
#include "hictkpy/reference.hpp"

namespace nb = nanobind;

namespace hictkpy {

HiCFileWriter::HiCFileWriter(std::string_view path, const nb::dict &chromosomes,
                             const std::vector<std::uint32_t> &resolutions,
                             std::string_view assembly, std::size_t n_threads,
                             std::size_t chunk_size, const std::filesystem::path &tmpdir,
                             std::uint32_t compression_lvl, bool skip_all_vs_all_matrix)
    : _tmpdir(tmpdir, true),
      _w(path, chromosome_dict_to_reference(chromosomes), resolutions, assembly, n_threads,
         chunk_size, _tmpdir(), compression_lvl, skip_all_vs_all_matrix) {}

HiCFileWriter::HiCFileWriter(std::string_view path, const nb::dict &chromosomes,
                             std::uint32_t resolution, std::string_view assembly,
                             std::size_t n_threads, std::size_t chunk_size,
                             const std::filesystem::path &tmpdir, std::uint32_t compression_lvl,
                             bool skip_all_vs_all_matrix)
    : HiCFileWriter(path, chromosomes, std::vector<std::uint32_t>{resolution}, assembly, n_threads,
                    chunk_size, tmpdir, compression_lvl, skip_all_vs_all_matrix) {}

void HiCFileWriter::serialize([[maybe_unused]] const std::string &log_lvl_str) {
  if (_finalized) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("finalize was already called on file \"{}\""), _w.path()));
  }
#ifndef _WIN32
  // TODO fixme
  // There is something very odd going on when trying to call most spdlog functions from within
  // Python bindings on recent versions of Windows.
  // Possibly related to https://github.com/gabime/spdlog/issues/3212
  const auto log_lvl = spdlog::level::from_str(log_lvl_str);
  const auto previous_lvl = spdlog::default_logger()->level();
  spdlog::default_logger()->set_level(log_lvl);
#endif
  _w.serialize();
  _finalized = true;
#ifndef _WIN32
  spdlog::default_logger()->set_level(previous_lvl);
#endif
}

std::string_view HiCFileWriter::path() const noexcept { return _w.path(); }

const std::vector<std::uint32_t> &HiCFileWriter::resolutions() const noexcept {
  return _w.resolutions();
}

const hictk::Reference &HiCFileWriter::chromosomes() const { return _w.chromosomes(); }

void HiCFileWriter::add_pixels(const nb::object &df) {
  const auto coo_format = nb::cast<bool>(df.attr("columns").attr("__contains__")("bin1_id"));
  const auto pixels =
      coo_format ? coo_df_to_thin_pixels<float>(df, false)
                 : bg2_df_to_thin_pixels<float>(_w.bins(_w.resolutions().front()), df, false);
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
  writer.def(nb::init<std::string_view, nb::dict, std::uint32_t, std::string_view, std::size_t,
                      std::size_t, const std::filesystem::path &, std::uint32_t, bool>(),
             nb::arg("path"), nb::arg("chromosomes"), nb::arg("resolution"),
             nb::arg("assembly") = "unknown", nb::arg("n_threads") = 1,
             nb::arg("chunk_size") = 10'000'000,
             nb::arg("tmpdir") = std::filesystem::temp_directory_path().string(),
             nb::arg("compression_lvl") = 10, nb::arg("skip_all_vs_all_matrix") = false,
             "Open a .hic file for writing.");

  writer.def(
      nb::init<std::string_view, nb::dict, const std::vector<std::uint32_t> &, std::string_view,
               std::size_t, std::size_t, const std::filesystem::path &, std::uint32_t, bool>(),
      nb::arg("path"), nb::arg("chromosomes"), nb::arg("resolutions"),
      nb::arg("assembly") = "unknown", nb::arg("n_threads") = 1, nb::arg("chunk_size") = 10'000'000,
      nb::arg("tmpdir") = std::filesystem::temp_directory_path().string(),
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
  writer.def("finalize", &hictkpy::HiCFileWriter::serialize, nb::arg("log_lvl") = "warn",
             "Write interactions to file.");
}

}  // namespace hictkpy
