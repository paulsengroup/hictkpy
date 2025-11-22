// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/multires_file.hpp"

#include <fmt/format.h>

#include <cstdint>
#include <exception>
#include <hictk/cooler/multires_cooler.hpp>
#include <hictk/cooler/validation.hpp>
#include <hictk/hic.hpp>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictkpy/file.hpp"
#include "hictkpy/nanobind.hpp"
#include "hictkpy/reference.hpp"
#include "hictkpy/to_numpy.hpp"

namespace nb = nanobind;

namespace hictkpy {

static std::string repr(const MultiResFile& mrf) {
  return fmt::format(FMT_STRING("MultiResFile({})"), mrf->path());
}

static MultiResFile& ctx_enter(MultiResFile& mrf) { return mrf; }

static void ctx_exit(MultiResFile& mrf, [[maybe_unused]] nb::handle exc_type,
                     [[maybe_unused]] nb::handle exc_value, [[maybe_unused]] nb::handle traceback) {
  mrf.try_close();
}

// macro used to reduce boilerplate required to call methods on File objects
#define HICTKPY_CALL_METHOD_CHECKED(method) [](const MultiResFile& x) { return x->method(); }

static File getitem(const MultiResFile& mrf, std::int64_t resolution) {
  constexpr auto max_resolution =
      static_cast<std::int64_t>(std::numeric_limits<std::uint32_t>::max());

  if (resolution < 0) {
    throw std::invalid_argument("resolution must be non-negative");
  }
  if (resolution > max_resolution) {
    throw std::invalid_argument(
        fmt::format(FMT_STRING("resolution must be less than {}"), max_resolution + 1));
  }

  return File{mrf->open(static_cast<std::uint32_t>(resolution))};
}

static auto get_chromosomes(const MultiResFile& mrf, bool include_ALL) {
  return get_chromosomes_from_object(*mrf, include_ALL);
}

static std::filesystem::path get_path(const MultiResFile& mrf) {
  return std::filesystem::path{mrf->path()};
}

static auto get_resolutions(const MultiResFile& mrf) {
  return make_owning_numpy<std::int64_t>(mrf->resolutions());
}

static nb::dict get_attrs(const hictk::hic::File& hf) {
  nb::dict py_attrs;

  py_attrs["format"] = "HIC";
  py_attrs["format-version"] = hf.version();
  py_attrs["assembly"] = hf.assembly();
  py_attrs["format-url"] = "https://github.com/aidenlab/hic-format";
  py_attrs["nchroms"] = hf.nchroms();

  for (const auto& [k, v] : hf.attributes()) {
    py_attrs[nb::cast(k)] = v;
  }

  return py_attrs;
}

static nb::dict get_attrs(const hictk::cooler::MultiResFile& mclr) {
  nb::dict py_attrs;

  py_attrs["format"] = mclr.attributes().format;
  py_attrs["format-version"] = mclr.attributes().format_version;
  py_attrs["format-url"] = "https://github.com/open2c/cooler";
  py_attrs["assembly"] =
      mclr.open(mclr.resolutions().front()).attributes().assembly.value_or("unknown");
  py_attrs["nchroms"] = mclr.chromosomes().size();

  return py_attrs;
}

static nb::dict get_attributes(const MultiResFile& mrf) {
  auto attrs = mrf->is_hic()
                   ? get_attrs(mrf->open(mrf->resolutions().front()).get<hictk::hic::File>())
                   : get_attrs(hictk::cooler::MultiResFile{mrf->path()});
  attrs["resolutions"] = get_resolutions(mrf);

  return attrs;
}

[[noreturn]] static void throw_closed_file_exc(std::string_view path) {
  throw std::runtime_error(fmt::format(
      FMT_STRING("caught an attempt to access file \"{}\", which has already been closed"), path));
}

MultiResFile::MultiResFile(const std::filesystem::path& path)
    : _fp(path.string()), _uri(path.string()) {}

hictk::MultiResFile* MultiResFile::operator->() { return &**this; }

const hictk::MultiResFile* MultiResFile::operator->() const { return &**this; }

hictk::MultiResFile& MultiResFile::operator*() {
  if (_fp.has_value()) {
    return *_fp;
  }
  throw_closed_file_exc(_uri);
}

const hictk::MultiResFile& MultiResFile::operator*() const {
  if (_fp.has_value()) {
    return *_fp;
  }
  throw_closed_file_exc(_uri);
}

void MultiResFile::close() { _fp.reset(); }

bool MultiResFile::try_close() noexcept {
  if (!_fp) {
    return true;
  }

  try {
    try {
      close();
      return true;
    } catch (const std::exception& e) {
      nanobind::module_::import_("warnings").attr("warn")(e.what());
    } catch (...) {
      nanobind::module_::import_("warnings")
          .attr("warn")(fmt::format(
              FMT_STRING("an error occurred while closing file \"{}\": unknown error"), _uri));
    }
  } catch (...) {  // NOLINT
  }

  return false;
}

bool MultiResFile::is_mcool(const std::filesystem::path& path) {
  return bool(hictk::cooler::utils::is_multires_file(path.string()));
}

void MultiResFile::bind(nb::module_& m) {
  auto mres_file = nb::class_<MultiResFile>(
      m, "MultiResFile", "Class representing a file handle to a .hic or .mcool file");
  mres_file.def(nb::init<const std::filesystem::path&>(), nb::arg("path"),
                "Open a multi-resolution Cooler file (.mcool) or .hic file.");

  mres_file.def("__repr__", &repr, nb::rv_policy::move);

  mres_file.def("__enter__", &ctx_enter, nb::rv_policy::reference_internal);

  mres_file.def("__exit__", &ctx_exit,
                // clang-format off
                nb::arg("exc_type") = nb::none(),
                nb::arg("exc_value") = nb::none(),
                nb::arg("traceback") = nb::none()
                // clang-format on
  );

  mres_file.def("path", &get_path, "Get the file path.", nb::rv_policy::move);
  mres_file.def("is_mcool", HICTKPY_CALL_METHOD_CHECKED(is_mcool),
                "Test whether the file is in .mcool format.");
  mres_file.def("is_hic", HICTKPY_CALL_METHOD_CHECKED(is_hic),
                "Test whether the file is in .hic format.");

  mres_file.def("close", &MultiResFile::close, "Manually close the file handle.");

  mres_file.def("chromosomes", &get_chromosomes, nb::arg("include_ALL") = false,
                "Get the chromosome sizes as a dictionary mapping names to sizes.",
                nb::rv_policy::take_ownership);
  mres_file.def("resolutions", &get_resolutions, "Get the list of available resolutions.",
                nb::rv_policy::take_ownership);
  mres_file.def("attributes", &get_attributes, "Get file attributes as a dictionary.",
                nb::rv_policy::take_ownership);
  mres_file.def("__getitem__", &getitem,
                "Open the Cooler or .hic file corresponding to the resolution given as input.",
                nb::rv_policy::move);
}

}  // namespace hictkpy
