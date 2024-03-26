// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <filesystem>
#include <hictk/bin_table.hpp>
#include <hictk/common.hpp>
#include <hictk/pixel.hpp>
#include <stdexcept>
#include <string_view>

#include "hictkpy/common.hpp"
#include "hictkpy/file_creation.hpp"

namespace nb = nanobind;

namespace hictkpy {

template <typename N>
std::vector<hictk::ThinPixel<N>> coo_df_to_thin_pixels(nanobind::object df, bool sorted) {
  using BufferT1 = nb::ndarray<nb::numpy, nb::shape<nb::any>, std::uint64_t>;
  using BufferT2 = nb::ndarray<nb::numpy, nb::shape<nb::any>, N>;

  auto bin1_ids_np = nb::cast<BufferT1>(df.attr("__getitem__")("bin1_id").attr("to_numpy")());
  auto bin2_ids_np = nb::cast<BufferT1>(df.attr("__getitem__")("bin2_id").attr("to_numpy")());
  auto counts_np = nb::cast<BufferT2>(df.attr("__getitem__")("count").attr("to_numpy")());

  const auto bin1_ids = bin1_ids_np.view();
  const auto bin2_ids = bin2_ids_np.view();
  const auto counts = counts_np.view();

  std::vector<hictk::ThinPixel<N>> buffer(bin1_ids_np.size());
  for (std::size_t i = 0; i < bin1_ids_np.size(); ++i) {
    buffer[i] = hictk::ThinPixel<N>{bin1_ids(i), bin2_ids(i), counts(i)};
  }

  if (sorted) {
    std::sort(buffer.begin(), buffer.end());
  }

  return buffer;
}

template <typename N>
std::vector<hictk::ThinPixel<N>> bg2_df_to_thin_pixels(const hictk::BinTable &bin_table,
                                                       nanobind::object df, bool sorted) {
  using BufferT1 = nb::ndarray<nb::numpy, nb::shape<nb::any>, std::uint32_t>;
  using BufferT2 = nb::ndarray<nb::numpy, nb::shape<nb::any>, N>;

  auto chrom1 = nb::cast<nb::list>(df.attr("__getitem__")("chrom1").attr("tolist")());
  auto start1_np = nb::cast<BufferT1>(df.attr("__getitem__")("start1").attr("to_numpy")());
  auto end1_np = nb::cast<BufferT1>(df.attr("__getitem__")("end1").attr("to_numpy")());
  auto chrom2 = nb::cast<nb::list>(df.attr("__getitem__")("chrom2").attr("tolist")());
  auto start2_np = nb::cast<BufferT1>(df.attr("__getitem__")("start2").attr("to_numpy")());
  auto end2_np = nb::cast<BufferT1>(df.attr("__getitem__")("end2").attr("to_numpy")());
  auto counts_np = nb::cast<BufferT2>(df.attr("__getitem__")("count").attr("to_numpy")());

  const auto start1 = start1_np.view();
  const auto end1 = end1_np.view();
  const auto start2 = start2_np.view();
  const auto end2 = end2_np.view();
  const auto counts = counts_np.view();

  const auto &reference = bin_table.chromosomes();
  std::vector<hictk::ThinPixel<N>> buffer(start1_np.size());
  for (std::size_t i = 0; i < start1_np.size(); ++i) {
    auto bin1 = bin_table.at(reference.at(nb::cast<nb::str>(chrom1[i]).c_str()), start1(i));
    auto bin2 = bin_table.at(reference.at(nb::cast<nb::str>(chrom2[i]).c_str()), start2(i));

    if (end1(i) - start1(i) > bin_table.resolution() ||
        end2(i) - start2(i) > bin_table.resolution()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Found an invalid pixel {} {} {} {} {} {} {}: pixel spans a "
                                 "distance greater than the bin size"),
                      nb::cast<nb::str>(chrom1[i]).c_str(), start1(i), end1(i),
                      nb::cast<nb::str>(chrom2[i]).c_str(), start2(i), end2(i), counts(i)));
    }

    buffer[i] =
        hictk::Pixel<N>{hictk::PixelCoordinates{std::move(bin1), std::move(bin2)}, counts(i)}
            .to_thin();
  }

  if (sorted) {
    std::sort(buffer.begin(), buffer.end());
  }

  return buffer;
}

HiCFileWriter::HiCFileWriter(std::string_view path, hictk::Reference chromosomes,
                             const std::vector<std::uint32_t> &resolutions,
                             std::string_view assembly, std::size_t n_threads,
                             std::size_t chunk_size, const std::filesystem::path &tmpdir,
                             std::uint32_t compression_lvl, bool skip_all_vs_all_matrix)
    : _w(path, std::move(chromosomes), resolutions, assembly, n_threads, chunk_size, tmpdir,
         compression_lvl, skip_all_vs_all_matrix) {}

void HiCFileWriter::serialize(const std::string &log_lvl_str) {
  if (_finalized) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("finalize was already called on file \"{}\""), _w.path()));
  }
  spdlog::level::level_enum log_lvl = spdlog::level::from_str(log_lvl_str);
  const auto previous_lvl = spdlog::default_logger()->level();
  spdlog::default_logger()->set_level(log_lvl);
  _w.serialize();
  _finalized = true;
  spdlog::default_logger()->set_level(previous_lvl);
}

std::string_view HiCFileWriter::path() const noexcept { return _w.path(); }

const std::vector<std::uint32_t> &HiCFileWriter::resolutions() const noexcept {
  return _w.resolutions();
}

const hictk::Reference &HiCFileWriter::chromosomes() const { return _w.chromosomes(); }

void HiCFileWriter::add_pixels(nanobind::object df) {
  const auto coo_format = nb::cast<bool>(df.attr("columns").attr("__contains__")("bin1_id"));
  const auto pixels =
      coo_format ? coo_df_to_thin_pixels<float>(df, false)
                 : bg2_df_to_thin_pixels<float>(_w.bins(_w.resolutions().front()), df, false);
  _w.add_pixels(_w.resolutions().front(), pixels.begin(), pixels.end());
}

void hic_file_writer_ctor(hictkpy::HiCFileWriter *fp, std::string_view path,
                          nanobind::dict chromosomes, const std::vector<std::uint32_t> &resolutions,
                          std::string_view assembly, std::size_t n_threads, std::size_t chunk_size,
                          std::string_view tmpdir, std::uint32_t compression_lvl,
                          bool skip_all_vs_all_matrix) {
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  std::for_each(chromosomes.begin(), chromosomes.end(), [&](const auto &kv) {
    chrom_names.push_back(nb::cast<std::string>(kv.first));
    chrom_sizes.push_back(nb::cast<std::uint32_t>(kv.second));
  });

  new (fp)
      HiCFileWriter{path,
                    hictk::Reference{chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()},
                    resolutions,
                    assembly,
                    n_threads,
                    chunk_size,
                    std::filesystem::path{tmpdir},
                    compression_lvl,
                    skip_all_vs_all_matrix};
}

void hic_file_writer_ctor_single_res(hictkpy::HiCFileWriter *fp, std::string_view path,
                                     nanobind::dict chromosomes, std::uint32_t resolution,
                                     std::string_view assembly, std::size_t n_threads,
                                     std::size_t chunk_size, std::string_view tmpdir,
                                     std::uint32_t compression_lvl, bool skip_all_vs_all_matrix) {
  return hic_file_writer_ctor(fp, path, chromosomes, std::vector<std::uint32_t>{resolution},
                              assembly, n_threads, chunk_size, tmpdir, compression_lvl,
                              skip_all_vs_all_matrix);
}

std::string hic_file_writer_repr(hictkpy::HiCFileWriter &w) {
  return fmt::format(FMT_STRING("HiCFileWriter({})"), w.path());
}

CoolFileWriter::CoolFileWriter(std::string_view path_, hictk::Reference chromosomes_,
                               std::uint32_t resolution_, std::string_view assembly,
                               const std::filesystem::path &tmpdir, std::uint32_t compression_lvl)
    : _path(std::string{path_}),
      _tmpdir(tmpdir / (std::string{path_} + ".tmp")),
      _w(create_file(path_, std::move(chromosomes_), resolution_, _tmpdir())) {
  if (std::filesystem::exists(_path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to create .cool file \"{}\": file already exists"), path()));
  }
}

std::string CoolFileWriter::path() const noexcept { return _w->path(); }
std::uint32_t CoolFileWriter::resolution() const noexcept { return _w->resolution(); }
const hictk::Reference &CoolFileWriter::chromosomes() const { return _w->chromosomes(); }

void CoolFileWriter::add_pixels(nanobind::object df) {
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

void CoolFileWriter::serialize(const std::string &log_lvl_str) {
  if (_finalized) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("finalize was already called on file \"{}\""), _w->path()));
  }
  spdlog::level::level_enum log_lvl = spdlog::level::from_str(log_lvl_str);
  const auto previous_lvl = spdlog::default_logger()->level();
  spdlog::default_logger()->set_level(log_lvl);
  std::visit(
      [&](const auto &num) {
        using N = hictk::remove_cvref_t<decltype(num)>;
        _w->aggregate<N>(_path, false, _compression_lvl);
      },
      _w->open("0").pixel_variant());

  _finalized = true;
  spdlog::default_logger()->set_level(previous_lvl);
  const std::string sclr_path{_w->path()};
  _w.reset();
  std::filesystem::remove(sclr_path);
}

hictk::cooler::SingleCellFile CoolFileWriter::create_file(std::string_view path,
                                                          hictk::Reference chromosomes,
                                                          std::uint32_t resolution,
                                                          const std::filesystem::path &tmpdir) {
  return hictk::cooler::SingleCellFile::create(tmpdir / std::filesystem::path{path}.filename(),
                                               std::move(chromosomes), resolution, false);
}

void cool_file_writer_ctor(hictkpy::CoolFileWriter *fp, std::string_view path,
                           nanobind::dict chromosomes, std::uint32_t resolution,
                           std::string_view assembly, std::string_view tmpdir,
                           std::uint32_t compression_lvl) {
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  std::for_each(chromosomes.begin(), chromosomes.end(), [&](const auto &kv) {
    chrom_names.push_back(nb::cast<std::string>(kv.first));
    chrom_sizes.push_back(nb::cast<std::uint32_t>(kv.second));
  });

  new (fp)
      CoolFileWriter{path,
                     hictk::Reference{chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()},
                     resolution,
                     assembly,
                     std::filesystem::path{tmpdir},
                     compression_lvl};
}

std::string cool_file_writer_repr(hictkpy::CoolFileWriter &w) {
  return fmt::format(FMT_STRING("CoolFileWriter({})"), w.path());
}

}  // namespace hictkpy
