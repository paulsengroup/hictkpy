// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>

#include "hictk/cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/multires_file.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/file.hpp"
#include "hictkpy/file_creation.hpp"
#include "hictkpy/multires_file.hpp"
#include "hictkpy/pixel_selector.hpp"
#include "hictkpy/singlecell_file.hpp"

namespace nb = nanobind;
namespace hictkpy {

template <typename N>
static void declare_thin_pixel_class(nb::module_ &m, const std::string &suffix) {
  const auto type_name = std::string{"ThinPixel"} + suffix;
  nb::class_<hictk::ThinPixel<N>>(m, type_name.c_str(), "Pixel in COO format.")
      .def_prop_ro("bin1_id", [](const hictk::ThinPixel<N> &tp) { return tp.bin1_id; })
      .def_prop_ro("bin2_id", [](const hictk::ThinPixel<N> &tp) { return tp.bin2_id; })
      .def_prop_ro("count", [](const hictk::ThinPixel<N> &tp) { return tp.count; })
      .def("__repr__",
           [](const hictk::ThinPixel<N> &tp) {
             return fmt::format(FMT_COMPILE("bin1_id={}; bin2_id={}; count={};"), tp.bin1_id,
                                tp.bin2_id, tp.count);
           })
      .def("__str__", [](const hictk::ThinPixel<N> &tp) {
        return fmt::format(FMT_COMPILE("{}\t{}\t{}"), tp.bin1_id, tp.bin2_id, tp.count);
      });
}

template <typename N>
static void declare_pixel_class(nb::module_ &m, const std::string &suffix) {
  const auto type_name = std::string{"Pixel"} + suffix;
  nb::class_<hictk::Pixel<N>>(m, type_name.c_str(), "Pixel in BG2 format.")
      .def_prop_ro("bin1_id", [](const hictk::Pixel<N> &p) { return p.coords.bin1.id(); })
      .def_prop_ro("bin2_id", [](const hictk::Pixel<N> &p) { return p.coords.bin2.id(); })
      .def_prop_ro("rel_bin1_id", [](const hictk::Pixel<N> &p) { return p.coords.bin1.rel_id(); })
      .def_prop_ro("rel_bin2_id", [](const hictk::Pixel<N> &p) { return p.coords.bin2.rel_id(); })
      .def_prop_ro("chrom1", [](const hictk::Pixel<N> &p) { return p.coords.bin1.chrom().name(); })
      .def_prop_ro("start1", [](const hictk::Pixel<N> &p) { return p.coords.bin1.start(); })
      .def_prop_ro("end1", [](const hictk::Pixel<N> &p) { return p.coords.bin1.end(); })
      .def_prop_ro("chrom2", [](const hictk::Pixel<N> &p) { return p.coords.bin2.chrom().name(); })
      .def_prop_ro("start2", [](const hictk::Pixel<N> &p) { return p.coords.bin2.start(); })
      .def_prop_ro("end2", [](const hictk::Pixel<N> &p) { return p.coords.bin2.end(); })
      .def_prop_ro("count", [](const hictk::Pixel<N> &p) { return p.count; })
      .def("__repr__",
           [](const hictk::Pixel<N> &p) {
             return fmt::format(
                 FMT_COMPILE("chrom1={}; start1={}; end1={}; chrom2={}; start2={}; end2={};"),
                 p.coords.bin1.chrom().name(), p.coords.bin1.start(), p.coords.bin1.end(),
                 p.coords.bin2.chrom().name(), p.coords.bin2.start(), p.coords.bin2.end(), p.count);
           })
      .def("__str__", [](const hictk::Pixel<N> &p) {
        return fmt::format(FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{}"), p.coords.bin1.chrom().name(),
                           p.coords.bin1.start(), p.coords.bin1.end(), p.coords.bin2.chrom().name(),
                           p.coords.bin2.start(), p.coords.bin2.end(), p.count);
      });
}

static void declare_pixel_selector_class(nb::module_ &m) {
  auto sel = nb::class_<PixelSelector>(
      m, "PixelSelector",
      "Class representing pixels overlapping with the given genomic intervals.");

  sel.def(nb::init<std::shared_ptr<const hictk::cooler::PixelSelector>, std::string_view, bool, bool>(),
          nb::arg("selector"), nb::arg("type"), nb::arg("join"), nb::arg("_mirror"));
  sel.def(nb::init<std::shared_ptr<const hictk::hic::PixelSelector>, std::string_view, bool, bool>(),
          nb::arg("selector"), nb::arg("type"), nb::arg("join"), nb::arg("_mirror"));
  sel.def(nb::init<std::shared_ptr<const hictk::hic::PixelSelectorAll>, std::string_view, bool, bool>(),
          nb::arg("selector"), nb::arg("type"), nb::arg("join"), nb::arg("_mirror"));

  sel.def("__repr__", &PixelSelector::repr);

  sel.def("coord1", &PixelSelector::get_coord1, "Get query coordinates for the first dimension.");
  sel.def("coord2", &PixelSelector::get_coord2, "Get query coordinates for the second dimension.");

  sel.def("__iter__", &PixelSelector::make_iterable, nb::keep_alive<0, 1>());

  sel.def("to_df", &PixelSelector::to_df, "Retrieve interactions as a pandas DataFrame.");
  sel.def("to_numpy", &PixelSelector::to_numpy, "Retrieve interactions as a numpy 2D matrix.");
  sel.def("to_coo", &PixelSelector::to_coo, "Retrieve interactions as a scipy.sparse.coo_matrix.");

  sel.def("nnz", &PixelSelector::nnz,
          "Get the number of non-zero entries for the current pixel selection.");
  sel.def("sum", &PixelSelector::sum,
          "Get the total number of interactions for the current pixel selection.");
}

static void declare_file_class(nb::module_ &m) {
  auto file = nb::class_<hictk::File>(m, "File",
                                      "Class representing a file handle to a .cool or .hic file.");

  file.def("__init__", &file::ctor, nb::arg("path"), nb::arg("resolution"),
           nb::arg("matrix_type") = "observed", nb::arg("matrix_unit") = "BP",
           "Construct a file object to a .hic, .cool or .mcool file given the file path and "
           "resolution.\n"
           "Resolution is ignored when opening single-resolution Cooler files.");

  file.def("__repr__", &file::repr);

  file.def("uri", &hictk::File::uri, "Return the file URI.");
  file.def("path", &hictk::File::path, "Return the file path.");

  file.def("is_hic", &hictk::File::is_hic, "Test whether file is in .hic format.");
  file.def("is_cooler", &hictk::File::is_cooler, "Test whether file is in .cool format.");

  file.def("chromosomes", &get_chromosomes_from_file<hictk::File>, nb::arg("include_all") = false,
           "Get chromosomes sizes as a dictionary mapping names to sizes.");
  file.def("bins", &get_bins_from_file<hictk::File>, "Get bins as a pandas DataFrame.");

  file.def("resolution", &hictk::File::resolution, "Get the bin size in bp.");
  file.def("nbins", &hictk::File::nbins, "Get the total number of bins.");
  file.def("nchroms", &hictk::File::nchroms, "Get the total number of chromosomes.");

  file.def("attributes", &file::attributes, "Get file attributes as a dictionary.");

  file.def("fetch", &file::fetch, nb::keep_alive<0, 1>(), nb::arg("range1") = "",
           nb::arg("range2") = "", nb::arg("normalization") = "NONE", nb::arg("count_type") = "int",
           nb::arg("join") = false, nb::arg("query_type") = "UCSC",
           "Fetch interactions overlapping a region of interest.");

  file.def("avail_normalizations", &file::avail_normalizations,
           "Get the list of available normalizations.");
  file.def("has_normalization", &hictk::File::has_normalization, nb::arg("normalization"),
           "Check whether a given normalization is available.");
}

static void declare_multires_file_class(nb::module_ &m) {
  auto mres_file = nb::class_<hictk::MultiResFile>(
      m, "MultiResFile", "Class representing a file handle to a .hic or .mcool file");
  mres_file.def("__init__", &multires_file::ctor, nb::arg("path"),
                "Open a multi-resolution Cooler file (.mcool).");

  mres_file.def("__repr__", &multires_file::repr);

  mres_file.def("path", &hictk::MultiResFile::path, "Get the file path.");
  mres_file.def("chromosomes", &get_chromosomes_from_file<hictk::MultiResFile>,
                nb::arg("include_all") = false,
                "Get chromosomes sizes as a dictionary mapping names to sizes.");
  mres_file.def("resolutions", &hictk::MultiResFile::resolutions,
                "Get the list of available resolutions.");
  mres_file.def("__getitem__", &hictk::MultiResFile::open,
                "Open the Cooler file corresponding to the resolution given as input.");
}

static void declare_singlecell_file_class(nb::module_ &m) {
  auto cooler = m.def_submodule("cooler");

  auto scell_file = nb::class_<hictk::cooler::SingleCellFile>(
      cooler, "SingleCellFile", "Class representing a file handle to a .scool file.");

  scell_file.def("__init__", &singlecell_file::ctor, nb::arg("path"),
                 "Open a single-cell Cooler file (.scool).");

  scell_file.def("__repr__", &singlecell_file::repr);

  scell_file.def("path", &hictk::cooler::SingleCellFile::path, "Get the file path.");
  scell_file.def("resolution", &hictk::cooler::SingleCellFile::resolution,
                 "Get the bin size in bp.");
  scell_file.def("chromosomes", &get_chromosomes_from_file<hictk::cooler::SingleCellFile>,
                 nb::arg("include_all") = false,
                 "Get chromosomes sizes as a dictionary mapping names to sizes.");
  scell_file.def("bins", &get_bins_from_file<hictk::cooler::SingleCellFile>,
                 "Get bins as a pandas DataFrame.");
  scell_file.def("attributes", &singlecell_file::get_attrs, "Get file attributes as a dictionary.");
  scell_file.def("cells", &singlecell_file::get_cells, "Get the list of available cells.");
  scell_file.def("__getitem__", &singlecell_file::getitem,
                 "Open the Cooler file corresponding to the cell ID given as input.");
}

static void declare_hic_file_writer_class(nb::module_ &m) {
  auto hic = m.def_submodule("hic");

  auto writer = nb::class_<hictkpy::HiCFileWriter>(
      hic, "FileWriter", "Class representing a file handle to create .hic files.");

  writer.def("__init__", &hic_file_writer_ctor_single_res, nb::arg("path"), nb::arg("chromosomes"),
             nb::arg("resolution"), nb::arg("assembly") = "unknown", nb::arg("n_threads") = 1,
             nb::arg("chunk_size") = 10'000'000,
             nb::arg("tmpdir") = std::filesystem::temp_directory_path().string(),
             nb::arg("compression_lvl") = 9, nb::arg("skip_all_vs_all_matrix") = false,
             "Open a .hic file for writing.");

  writer.def("__init__", &hic_file_writer_ctor, nb::arg("path"), nb::arg("chromosomes"),
             nb::arg("resolutions"), nb::arg("assembly") = "unknown", nb::arg("n_threads") = 1,
             nb::arg("chunk_size") = 10'000'000,
             nb::arg("tmpdir") = std::filesystem::temp_directory_path().string(),
             nb::arg("compression_lvl") = 9, nb::arg("skip_all_vs_all_matrix") = false,
             "Open a .hic file for writing.");

  writer.def("__repr__", &hic_file_writer_repr);

  writer.def("path", &hictkpy::HiCFileWriter::path, "Get the file path.");
  writer.def("resolutions", &hictkpy::HiCFileWriter::resolutions,
             "Get the list of resolutions in bp.");
  writer.def("chromosomes", &get_chromosomes_from_file<hictkpy::HiCFileWriter>,
             nb::arg("include_all") = false,
             "Get chromosomes sizes as a dictionary mapping names to sizes.");

  writer.def("add_pixels", &hictkpy::HiCFileWriter::add_pixels, nb::arg("pixels"),
             "Add pixels from a pandas DataFrame containing pixels in COO or BG2 format (i.e. "
             "either with columns=[bin1_id, bin2_id, count] or with columns=[chrom1, start1, end1, "
             "chrom2, start2, end2, count].");
  writer.def("finalize", &hictkpy::HiCFileWriter::serialize, nb::arg("log_lvl") = "warn",
             "Write interactions to file.");
}

static void declare_cooler_file_writer_class(nb::module_ &m) {
  auto cooler = m.def_submodule("cooler");

  auto writer = nb::class_<hictkpy::CoolFileWriter>(
      cooler, "FileWriter", "Class representing a file handle to create .cool files.");

  writer.def("__init__", &cool_file_writer_ctor, nb::arg("path"), nb::arg("chromosomes"),
             nb::arg("resolution"), nb::arg("assembly") = "unknown",
             nb::arg("tmpdir") = std::filesystem::temp_directory_path().string(),
             nb::arg("compression_lvl") = 6, "Open a .cool file for writing.");

  writer.def("__repr__", &cool_file_writer_repr);

  writer.def("path", &hictkpy::CoolFileWriter::path, "Get the file path.");
  writer.def("resolutions", &hictkpy::CoolFileWriter::resolution, "Get the resolution in bp.");
  writer.def("chromosomes", &get_chromosomes_from_file<hictkpy::CoolFileWriter>,
             nb::arg("include_all") = false,
             "Get chromosomes sizes as a dictionary mapping names to sizes.");

  writer.def("add_pixels", &hictkpy::CoolFileWriter::add_pixels, nb::arg("pixels"),
             "Add pixels from a pandas DataFrame containing pixels in COO or BG2 format (i.e. "
             "either with columns=[bin1_id, bin2_id, count] or with columns=[chrom1, start1, end1, "
             "chrom2, start2, end2, count].");
  writer.def("finalize", &hictkpy::CoolFileWriter::serialize, nb::arg("log_lvl") = "warn",
             "Write interactions to file.");
}

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_hictkpy, m) {
  [[maybe_unused]] auto np = nb::module_::import_("numpy");
  [[maybe_unused]] auto pd = nb::module_::import_("pandas");
  [[maybe_unused]] auto ss = nb::module_::import_("scipy.sparse");

  m.attr("__hictk_version__") = hictk::config::version::str();

  m.doc() = "Blazing fast toolkit to work with .hic and .cool files.";

  m.def("is_cooler", &file::is_cooler, nb::arg("path"),
        "Test whether path points to a cooler file.");
  m.def("is_mcool_file", &file::is_mcool_file, nb::arg("path"),
        "Test whether path points to a .mcool file.");
  m.def("is_scool_file", &file::is_scool_file, nb::arg("path"),
        "Test whether path points to a .scool file.");
  m.def("is_hic", &file::is_hic, nb::arg("path"), "Test whether path points to a .hic file.");

  declare_thin_pixel_class<std::int32_t>(m, "Int");
  declare_thin_pixel_class<double>(m, "FP");
  declare_pixel_class<std::int32_t>(m, "Int");
  declare_pixel_class<double>(m, "FP");

  declare_pixel_selector_class(m);

  declare_file_class(m);

  declare_multires_file_class(m);
  declare_singlecell_file_class(m);

  declare_hic_file_writer_class(m);
  declare_cooler_file_writer_class(m);
}

}  // namespace hictkpy
