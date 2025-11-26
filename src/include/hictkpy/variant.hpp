// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <variant>

namespace hictkpy {

// clang-format off
using NumericDtype =
    std::variant<
        std::uint8_t,
        std::uint16_t,
        std::uint32_t,
        std::uint64_t,
        std::int8_t,
        std::int16_t,
        std::int32_t,
        std::int64_t,
        float,
        double
    >;
// clang-format on

// clang-format off
using Dtype =
    std::variant<
        std::monostate,
        std::uint8_t,
        std::uint16_t,
        std::uint32_t,
        std::uint64_t,
        std::int8_t,
        std::int16_t,
        std::int32_t,
        std::int64_t,
        float,
        double,
        std::string
    >;
// clang-format on

}  // namespace hictkpy
