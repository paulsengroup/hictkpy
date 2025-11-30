// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/nanobind.hpp"

#include <Python.h>
#include <fmt/format.h>

#include <cassert>
#include <exception>
#include <hictk/numeric_utils.hpp>
#include <string>
#include <utility>
#include <vector>

#include "hictkpy/locking.hpp"

namespace nb = nanobind;

namespace hictkpy {

nb::module_ import_module_checked(const std::string& module_name) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  try {
    return nb::module_::import_(module_name.c_str());
  } catch (nb::python_error& e) {
    // NOLINTNEXTLINE(*-vararg)
    nb::raise_from(e, PyExc_ModuleNotFoundError,
                   "To enable %s support, please install %s with: pip install 'hictkpy[%s]'\n"
                   "Alternatively, you can install hictkpy with all its dependencies by "
                   "running: pip install 'hictkpy[all]'",
                   module_name.c_str(), module_name.c_str(), module_name.c_str());
  }
}

[[noreturn]] static void throw_import_error(const char* msg) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  throw nb::import_error(msg);
}

template <typename... T>
[[noreturn]] static void throw_import_error(fmt::format_string<T...> fmt, T&&... args) {
  const auto msg = fmt::format(fmt, std::forward<T>(args)...);
  throw_import_error(msg.c_str());
}

[[nodiscard]] static auto get_pyarrow_version() {
  struct PyArrowVersion {
    int major{};
    int minor{};
    int patch{};
  };

  [[maybe_unused]] const GilScopedAcquire gil{true};
  auto metadata = nb::module_::import_("importlib.metadata");

  const auto version =
      nb::cast<std::vector<std::string>>(metadata.attr("version")("pyarrow").attr("split")("."));
  if (version.size() < 3) {
    throw_import_error(
        "Unable to detect pyarrow version: assuming pyarrow's version is not compatible.\n"
        "Please install a compatible version of pyarrow with: pip install 'hictkpy[pyarrow]'");
  }

  return PyArrowVersion{hictk::internal::parse_numeric_or_throw<int>(version[0]),
                        hictk::internal::parse_numeric_or_throw<int>(version[1]),
                        hictk::internal::parse_numeric_or_throw<int>(version[2])};
}

template <typename Version>
[[nodiscard]] static bool check_version_is_new_enough(const Version& version_found,
                                                      int min_version_major, int min_version_minor,
                                                      int min_version_patch) {
  auto version_ok = version_found.major >= min_version_major;
  version_ok |=
      version_found.major == min_version_major && version_found.minor >= min_version_minor;
  version_ok |= version_found.major == min_version_major &&
                version_found.minor == min_version_minor &&
                version_found.patch >= min_version_patch;

  if (!version_ok) {
    throw_import_error(
        FMT_STRING("pyarrow {}.{}.{} is too old to be used with hictkpy.\n"
                   "Please install a compatible version with: pip install 'hictkpy[pyarrow]'"),
        version_found.major, version_found.minor, version_found.patch);
  }

  return true;
}

nb::module_ import_pyarrow_checked(int min_version_major, int min_version_minor,
                                   int min_version_patch) {
  assert(min_version_major >= 0);
  assert(min_version_minor >= 0);
  assert(min_version_patch >= 0);

  static bool version_ok{false};

  if (version_ok) {
    return import_module_checked("pyarrow");
  }

  try {
    version_ok = check_version_is_new_enough(get_pyarrow_version(), min_version_major,
                                             min_version_minor, min_version_patch);
    assert(version_ok);
    return import_module_checked("pyarrow");
  } catch (const nb::builtin_exception&) {
    throw;
  } catch (const std::exception& e) {
    throw_import_error(
        FMT_STRING(
            "Unable to parse pyarrow version: {}.\n"
            "Assuming pyarrow's version is not compatible.\n"
            "Please install a compatible version of pyarrow with: pip install 'hictkpy[pyarrow]'"),
        e.what());
  } catch (...) {
    throw_import_error(
        "Unable to parse pyarrow version.\n"
        "Assuming pyarrow's version is not compatible.\n"
        "Please install a compatible version of pyarrow with: pip install 'hictkpy[pyarrow]'");
  }
}

void check_module_is_importable(const std::string& module_name) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  std::ignore = import_module_checked(module_name);
}

void check_pyarrow_is_importable(int min_version_major, int min_version_minor,
                                 int min_version_patch) {
  [[maybe_unused]] const GilScopedAcquire gil{true};
  // NOLINTNEXTLINE(clang-analyzer-core.NullDereference)
  std::ignore = import_pyarrow_checked(min_version_major, min_version_minor, min_version_patch);
}

}  // namespace hictkpy
