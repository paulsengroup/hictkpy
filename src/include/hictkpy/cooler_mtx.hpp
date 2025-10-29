// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <mutex>

namespace hictkpy {

class CoolerGlobalLock {
  template <typename Mtx>
  class UniqueLock {
    std::unique_lock<Mtx> _lck{};

   public:
    UniqueLock() = default;
    UniqueLock(Mtx& mtx) : _lck(mtx) {  // NOLINT(*-explicit-conversions)
      SPDLOG_TRACE(FMT_STRING("CoolerGlobalLock({}): locked!"), fmt::ptr(&mtx));
    }

    UniqueLock(const UniqueLock&) = delete;
    UniqueLock(UniqueLock&&) noexcept = default;
    ~UniqueLock() noexcept {
      if (const auto* mtx = _lck.mutex(); !!mtx) {
        SPDLOG_TRACE(FMT_STRING("CoolerGlobalLock({}): unlocking..."), fmt::ptr(mtx));
      }
    }

    UniqueLock& operator=(const UniqueLock&) = delete;
    UniqueLock& operator=(UniqueLock&&) noexcept = default;

    [[nodiscard]] std::unique_lock<Mtx>& get() noexcept { return _lck; }
    [[nodiscard]] const std::unique_lock<Mtx>& get() const noexcept { return _lck; }
  };

  CoolerGlobalLock() = default;

 public:
  CoolerGlobalLock(const CoolerGlobalLock&) = delete;
  CoolerGlobalLock(CoolerGlobalLock&&) noexcept = delete;
  ~CoolerGlobalLock() noexcept = default;

  CoolerGlobalLock& operator=(const CoolerGlobalLock&) = delete;
  CoolerGlobalLock& operator=(CoolerGlobalLock&&) noexcept = delete;

  [[nodiscard]] static UniqueLock<std::recursive_mutex> lock() {
    static std::recursive_mutex mtx;
    SPDLOG_TRACE(FMT_STRING("CoolerGlobalLock({}): locking..."), fmt::ptr(&mtx));
    return {mtx};
  }
};

#define HICTKPY_LOCK_COOLER_MTX_SCOPED \
  [[maybe_unused]] const auto cooler_lock = CoolerGlobalLock::lock();

}  // namespace hictkpy
