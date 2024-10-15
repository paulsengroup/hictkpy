// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// Source: https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/
// GCC to MSVC codes: https://github.com/srz-zumix/awesome-cpp-warning

// clang-format off

// NOLINTBEGIN(cppcoreguidelines-macro-usage)

// Defines for MSVC
#ifdef _MSC_VER
    #define HICTKPY_DISABLE_WARNING_PUSH                      __pragma(warning(push))
    #define HICTKPY_DISABLE_WARNING_POP                       __pragma(warning(pop))
    #define HICTKPY_DISABLE_WARNING(warningNumber)            __pragma(warning(disable : warningNumber))

    #define HICTKPY_DISABLE_WARNING_CAST_ALIGN                HICTKPY_DISABLE_WARNING(4739)
    #define HICTKPY_DISABLE_WARNING_CXX98_COMPAT
    #define HICTKPY_DISABLE_WARNING_OLD_STYLE_CAST            HICTKPY_DISABLE_WARNING(4303)
    #define HICTKPY_DISABLE_WARNING_PEDANTIC                  HICTKPY_DISABLE_WARNING(4200) \
                                                              HICTKPY_DISABLE_WARNING(4201)
    #define HICTKPY_DISABLE_WARNING_SHADOW                    HICTKPY_DISABLE_WARNING(4456) \
                                                              HICTKPY_DISABLE_WARNING(4457) \
                                                              HICTKPY_DISABLE_WARNING(4458) \
                                                              HICTKPY_DISABLE_WARNING(4459)
    #define HICTKPY_DISABLE_WARNING_SIGN_CONVERSION           HICTKPY_DISABLE_WARNING(4308) \
                                                              HICTKPY_DISABLE_WARNING(4245) \
                                                              HICTKPY_DISABLE_WARNING(4365)
    #define HICTKPY_DISABLE_WARNING_USELESS_CAST

#endif

// Defines for GCC and Clang
#if defined(__GNUC__) || defined(__clang__)
    #define HICTKPY_DO_PRAGMA(X)                              _Pragma(#X)
    #define HICTKPY_DISABLE_WARNING_PUSH                      HICTKPY_DO_PRAGMA(GCC diagnostic push)
    #define HICTKPY_DISABLE_WARNING_POP                       HICTKPY_DO_PRAGMA(GCC diagnostic pop)
    #define HICTKPY_DISABLE_WARNING(warningName)              HICTKPY_DO_PRAGMA(GCC diagnostic ignored warningName)

    #define HICTKPY_DISABLE_WARNING_OLD_STYLE_CAST            HICTKPY_DISABLE_WARNING("-Wold-style-cast")
    #define HICTKPY_DISABLE_WARNING_PEDANTIC                  HICTKPY_DISABLE_WARNING("-Wpedantic")
    #define HICTKPY_DISABLE_WARNING_SHADOW                    HICTKPY_DISABLE_WARNING("-Wshadow")
    #define HICTKPY_DISABLE_WARNING_SIGN_CONVERSION           HICTKPY_DISABLE_WARNING("-Wsign-conversion")

#endif

// Defines for GCC only
#if defined(__GNUC__) && !defined(__clang__)
    #define HICTKPY_DISABLE_WARNING_CAST_ALIGN
    #define HICTKPY_DISABLE_WARNING_CXX98_COMPAT
    #define HICTKPY_DISABLE_WARNING_USELESS_CAST              HICTKPY_DISABLE_WARNING("-Wuseless-cast")
#endif

// Defines for Clang only
#ifdef __clang__
    #define HICTKPY_DISABLE_WARNING_CAST_ALIGN                HICTKPY_DISABLE_WARNING("-Wcast-align")
    #define HICTKPY_DISABLE_WARNING_CXX98_COMPAT              HICTKPY_DISABLE_WARNING("-Wc++98-compat-extra-semi")
    #define HICTKPY_DISABLE_WARNING_USELESS_CAST
#endif

// Defines for unknown/unsupported compilers
#if !defined(_MSC_VER) && !defined(__GNUC__) && !defined(__clang__)
    #define HICTKPY_DISABLE_WARNING
    #define HICTKPY_DISABLE_WARNING_PUSH
    #define HICTKPY_DISABLE_WARNING_POP

    #define HICTKPY_DISABLE_WARNING_CAST_ALIGN
    #define HICTKPY_DISABLE_WARNING_CXX98_COMPAT
    #define HICTKPY_DISABLE_WARNING_OLD_STYLE_CAST
    #define HICTKPY_DISABLE_WARNING_PEDANTIC
    #define HICTKPY_DISABLE_WARNING_SHADOW
    #define HICTKPY_DISABLE_WARNING_SIGN_CONVERSION
    #define HICTKPY_DISABLE_WARNING_USELESS_CAST

#endif

// NOLINTEND(cppcoreguidelines-macro-usage)

// clang-format on
