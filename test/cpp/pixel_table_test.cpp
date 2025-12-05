// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstdint>
#include <limits>

#include "hictkpy/pixel_table_helpers.hpp"

namespace hictkpy::test {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("hictkpy::PixelCounterCaster", "[table][short]") {
  using U8 = std::uint8_t;
  using U16 = std::uint16_t;
  using I8 = std::int8_t;
  using I16 = std::int16_t;
  SECTION("FP to FP") {
    CHECK(safe_numeric_cast<double>(float{}) == 0);
    CHECK(safe_numeric_cast<float>(double{}) == 0);
    CHECK(std::isnan(safe_numeric_cast<float>(std::numeric_limits<double>::quiet_NaN())));
    CHECK(std::isnan(safe_numeric_cast<double>(std::numeric_limits<float>::quiet_NaN())));
  }

  SECTION("FP to int") {
    CHECK(safe_numeric_cast<I8>(0.0) == 0);
    CHECK(safe_numeric_cast<I8>(0.1) == 0);
    CHECK(safe_numeric_cast<I8>(1.1) == 1);
    CHECK(safe_numeric_cast<I8>(-1.1) == -1);
    CHECK_THROWS(safe_numeric_cast<I8>(128));
    CHECK_THROWS(safe_numeric_cast<I8>(-129));
    CHECK_THROWS(safe_numeric_cast<I8>(std::numeric_limits<double>::quiet_NaN()));
  }

  SECTION("FP to uint") {
    CHECK(safe_numeric_cast<U8>(0.0) == 0);
    CHECK(safe_numeric_cast<U8>(0.1) == 0);
    CHECK(safe_numeric_cast<U8>(1.1) == 1);
    CHECK_THROWS(safe_numeric_cast<U8>(-1.0));
    CHECK_THROWS(safe_numeric_cast<U8>(256.0));
    CHECK_THROWS(safe_numeric_cast<U8>(std::numeric_limits<double>::quiet_NaN()));
  }

  SECTION("uint to uint") {
    CHECK(safe_numeric_cast<U8>(U8{1}) == 1);
    CHECK(safe_numeric_cast<U8>(U16{1}) == 1);
    CHECK_THROWS(safe_numeric_cast<U8>(U16{256}));
  }

  SECTION("uint to int") {
    CHECK(safe_numeric_cast<I8>(U8{1}) == 1);
    CHECK(safe_numeric_cast<I8>(U16{1}) == 1);
    CHECK_THROWS(safe_numeric_cast<I8>(U8{128}));
  }

  SECTION("uint to FP") { CHECK(safe_numeric_cast<double>(U8{0}) == 0); }

  SECTION("int to int") {
    CHECK(safe_numeric_cast<I8>(I8{1}) == 1);
    CHECK(safe_numeric_cast<I8>(I16{1}) == 1);
    CHECK(safe_numeric_cast<I8>(I16{-1}) == -1);
    CHECK_THROWS(safe_numeric_cast<I8>(I16{128}));
    CHECK_THROWS(safe_numeric_cast<I8>(I16{-129}));
  }

  SECTION("int to uint") {
    CHECK(safe_numeric_cast<U8>(I8{1}) == 1);
    CHECK_THROWS(safe_numeric_cast<U8>(I8{-1}));
    CHECK_THROWS(safe_numeric_cast<U8>(I16{256}));
  }

  SECTION("int to FP") { CHECK(safe_numeric_cast<double>(I8{0}) == 0); }

  SECTION("cast w/ name") {
    CHECK_THROWS_WITH(safe_numeric_cast<U8>("foo", I16{256}),
                      Catch::Matchers::Equals("unable to safely convert foo=256 (int16) to uint8"));
  }
}
// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictkpy::test
