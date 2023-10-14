// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include "hictk/cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictkpy/common.hpp"
#include "hictkpy/file.hpp"
#include "hictkpy/multires_file.hpp"
#include "hictkpy/pixel_selector.hpp"
#include "hictkpy/singlecell_file.hpp"

namespace py = pybind11;
namespace hictkpy {

template <typename N>
static void declare_thin_pixel_class(pybind11::module_ &m, const std::string &suffix) {
  const auto type_name = std::string{"ThinPixel"} + suffix;
  py::class_<hictk::ThinPixel<N>>(m, type_name.c_str())
      .def_property_readonly("bin1_id", [](const hictk::ThinPixel<N> &tp) { return tp.bin1_id; })
      .def_property_readonly("bin2_id", [](const hictk::ThinPixel<N> &tp) { return tp.bin2_id; })
      .def_property_readonly("count", [](const hictk::ThinPixel<N> &tp) { return tp.count; })
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
static void declare_pixel_class(pybind11::module_ &m, const std::string &suffix) {
  const auto type_name = std::string{"Pixel"} + suffix;
  py::class_<hictk::Pixel<N>>(m, type_name.c_str())
      .def_property_readonly("bin1_id", [](const hictk::Pixel<N> &p) { return p.coords.bin1.id(); })
      .def_property_readonly("bin2_id", [](const hictk::Pixel<N> &p) { return p.coords.bin2.id(); })
      .def_property_readonly("rel_bin1_id",
                             [](const hictk::Pixel<N> &p) { return p.coords.bin1.rel_id(); })
      .def_property_readonly("rel_bin2_id",
                             [](const hictk::Pixel<N> &p) { return p.coords.bin2.rel_id(); })
      .def_property_readonly("chrom1",
                             [](const hictk::Pixel<N> &p) { return p.coords.bin1.chrom().name(); })
      .def_property_readonly("start1",
                             [](const hictk::Pixel<N> &p) { return p.coords.bin1.start(); })
      .def_property_readonly("end1", [](const hictk::Pixel<N> &p) { return p.coords.bin1.end(); })
      .def_property_readonly("chrom2",
                             [](const hictk::Pixel<N> &p) { return p.coords.bin2.chrom().name(); })
      .def_property_readonly("start2",
                             [](const hictk::Pixel<N> &p) { return p.coords.bin2.start(); })
      .def_property_readonly("end2", [](const hictk::Pixel<N> &p) { return p.coords.bin2.end(); })
      .def_property_readonly("count", [](const hictk::Pixel<N> &p) { return p.count; })
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

static void declare_pixel_selector_class(pybind11::module_ &m) {
  auto sel =
      py::class_<PixelSelector>(m, "PixelSelector")
          .def(py::init<std::shared_ptr<const hictk::cooler::PixelSelector>, std::string_view,
                        bool>(),
               py::arg("selector"), py::arg("type"), py::arg("join"))
          .def(py::init<std::shared_ptr<const hictk::hic::PixelSelector>, std::string_view, bool>(),
               py::arg("selector"), py::arg("type"), py::arg("join"))
          .def(py::init<std::shared_ptr<const hictk::hic::PixelSelectorAll>, std::string_view,
                        bool>(),
               py::arg("selector"), py::arg("type"), py::arg("join"));

  sel.def("__repr__", &PixelSelector::repr);

  sel.def("coord1", &PixelSelector::get_coord1, "Get query coordinates for the first dimension.");
  sel.def("coord2", &PixelSelector::get_coord2, "Get query coordinates for the second dimension.");

  sel.def("__iter__", &PixelSelector::make_iterable, py::keep_alive<0, 1>());

  sel.def("to_df", &PixelSelector::to_df, "Retrieve interactions as a pandas DataFrame.");
  sel.def("to_numpy", &PixelSelector::to_numpy, "Retrieve interactions as a numpy 2D matrix.");
  sel.def("to_coo", &PixelSelector::to_coo, "Retrieve interactions as a scipy.sparse.coo_matrix.");

  sel.def("nnz", &PixelSelector::nnz,
          "Get the number of non-zero entries for the current pixel selection.");
  sel.def("sum", &PixelSelector::sum,
          "Get the total number of interactions for the current pixel selection.");
}

static void declare_file_class(pybind11::module_ &m) {
  auto file = py::class_<hictk::File>(m, "File").def(
      py::init(&file::ctor), py::arg("path"), py::arg("resolution"),
      py::arg("matrix_type") = "observed", py::arg("matrix_unit") = "BP",
      "Construct a file object to a .hic, .cool or .mcool file given the file path and "
      "resolution.\n"
      "Resolution is ignored when opening single-resolution Cooler files.");

  file.def("__repr__", &file::repr);

  file.def("uri", &hictk::File::uri, "Return the file URI.");
  file.def("path", &hictk::File::path, "Return the file path.");

  file.def("is_hic", &hictk::File::is_hic, "Test whether file is in .hic format.");
  file.def("is_cooler", &hictk::File::is_cooler, "Test whether file is in .cool format.");

  file.def("chromosomes", &get_chromosomes_from_file<hictk::File>, py::arg("include_all") = false,
           "Get chromosomes sizes as a dictionary mapping names to sizes.");
  file.def("bins", &get_bins_from_file<hictk::File>, "Get bins as a pandas DataFrame.");

  file.def("bin_size", &hictk::File::bin_size, "Get the bin size in bp.");
  file.def("nbins", &hictk::File::nbins, "Get the total number of bins.");
  file.def("nchroms", &hictk::File::nchroms, "Get the total number of chromosomes.");

  file.def("attributes", &file::attributes, "Get file attributes as a dictionary.");

  file.def("fetch", &file::fetch, py::keep_alive<0, 1>(), py::arg("range1") = "",
           py::arg("range2") = "", py::arg("normalization") = "NONE", py::arg("count_type") = "int",
           py::arg("join") = false, py::arg("query_type") = "UCSC",
           "Fetch interactions overlapping a region of interest.");

  file.def("avail_normalizations", &file::avail_normalizations,
           "Get the list of available normalizations.");
  file.def("has_normalization", &hictk::File::has_normalization, py::arg("normalization"),
           "Check whether a given normalization is available.");
}

static void declare_multires_file_class(pybind11::module_ &m) {
  auto cooler = m.def_submodule("cooler");

  auto mres_file = py::class_<hictk::cooler::MultiResFile>(cooler, "MultiResFile")
                       .def(py::init(&multires_file::ctor), py::arg("path"),
                            "Open a multi-resolution Cooler file (.mcool).");

  mres_file.def("__repr__", &multires_file::repr);

  mres_file.def("path", &hictk::cooler::MultiResFile::path, "Get the file path.");
  mres_file.def("chromosomes", &get_chromosomes_from_file<hictk::cooler::MultiResFile>,
                py::arg("include_all") = false,
                "Get chromosomes sizes as a dictionary mapping names to sizes.");
  mres_file.def("attributes", &multires_file::get_attrs, "Get file attributes as a dictionary.");
  mres_file.def("resolutions", &hictk::cooler::MultiResFile::resolutions,
                "Get the list of available resolutions.");
  mres_file.def("__getitem__", &multires_file::getitem,
                "Open the Cooler file corresponding to the resolution given as input.");
}

static void declare_singlecell_file_class(pybind11::module_ &m) {
  auto cooler = m.def_submodule("cooler");

  auto scell_file = py::class_<hictk::cooler::SingleCellFile>(cooler, "SingleCellFile")
                        .def(py::init(&singlecell_file::ctor), py::arg("path"),
                             "Open a single-cell Cooler file (.scool).");

  scell_file.def("__repr__", &singlecell_file::repr);

  scell_file.def("path", &hictk::cooler::SingleCellFile::path, "Get the file path.");
  scell_file.def("bin_size", &hictk::cooler::SingleCellFile::bin_size, "Get the bin size in bp.");
  scell_file.def("chromosomes", &get_chromosomes_from_file<hictk::cooler::SingleCellFile>,
                 py::arg("include_all") = false,
                 "Get chromosomes sizes as a dictionary mapping names to sizes.");
  scell_file.def("bins", &get_bins_from_file<hictk::cooler::SingleCellFile>,
                 "Get bins as a pandas DataFrame.");
  scell_file.def("attributes", &singlecell_file::get_attrs, "Get file attributes as a dictionary.");
  scell_file.def("cells", &singlecell_file::get_cells, "Get the list of available cells.");
  scell_file.def("__getitem__", &singlecell_file::getitem,
                 "Open the Cooler file corresponding to the cell ID given as input.");
}

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(hictkpy, m) {
  [[maybe_unused]] auto np = py::module::import("numpy");
  [[maybe_unused]] auto pd = py::module::import("pandas");
  [[maybe_unused]] auto ss = py::module::import("scipy.sparse");
  m.attr("__hictk_version__") = hictk::config::version::str();

  m.doc() = "Blazing fast toolkit to work with .hic and .cool files.";

  m.def("is_cooler", &file::is_cooler, py::arg("path"),
        "Test whether path points to a cooler file.");
  m.def("is_hic", &file::is_hic, py::arg("path"), "Test whether path points to a .hic file.");

  declare_thin_pixel_class<std::int32_t>(m, "Int");
  declare_thin_pixel_class<double>(m, "FP");
  declare_pixel_class<std::int32_t>(m, "Int");
  declare_pixel_class<double>(m, "FP");

  declare_pixel_selector_class(m);
  declare_file_class(m);
  declare_multires_file_class(m);
  declare_singlecell_file_class(m);
}

}  // namespace hictkpy
