// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <memory>
#include <mutex>

namespace hictkpy {

class CoolerGlobalLock {
  CoolerGlobalLock() = default;

 public:
  CoolerGlobalLock(const CoolerGlobalLock&) = delete;
  CoolerGlobalLock(CoolerGlobalLock&&) noexcept = delete;
  ~CoolerGlobalLock() noexcept = default;

  CoolerGlobalLock& operator=(const CoolerGlobalLock&) = delete;
  CoolerGlobalLock& operator=(CoolerGlobalLock&&) noexcept = delete;

  class UniqueLock {
   public:
    using Mutex = std::recursive_mutex;

   private:
    std::unique_lock<Mutex> _lck{};

   public:
    UniqueLock() = default;
    UniqueLock(Mutex& mtx) : _lck(mtx) {  // NOLINT(*-explicit-conversions)
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

    [[nodiscard]] std::unique_lock<Mutex>& get() noexcept { return _lck; }
    [[nodiscard]] const std::unique_lock<Mutex>& get() const noexcept { return _lck; }
  };

  [[nodiscard]] static UniqueLock::Mutex& mtx() {
    static UniqueLock::Mutex mtx_;
    return mtx_;
  }

  [[nodiscard]] static UniqueLock lock() {
    SPDLOG_TRACE(FMT_STRING("CoolerGlobalLock({}): locking..."), fmt::ptr(&mtx()));
    return {mtx()};
  }
};

#define HICTKPY_LOCK_COOLER_MTX_SCOPED \
  [[maybe_unused]] const auto cooler_lock = CoolerGlobalLock::lock();

class FileLock {
  using GlobalCoolerMutex = CoolerGlobalLock::UniqueLock::Mutex;
  std::shared_ptr<std::recursive_mutex> _mtx1{};
  std::unique_lock<std::recursive_mutex> _lck1{};
  std::unique_lock<GlobalCoolerMutex> _lck2{};

 public:
  FileLock() = default;
  FileLock(const FileLock&) = delete;
  FileLock(FileLock&&) noexcept = default;
  ~FileLock() noexcept = default;

  FileLock& operator=(const FileLock&) = delete;
  FileLock& operator=(FileLock&&) noexcept = default;

  explicit FileLock(std::shared_ptr<std::recursive_mutex> mtx, bool acquire = true);
  explicit FileLock(GlobalCoolerMutex& mtx, bool acquire = true);
  FileLock(std::shared_ptr<std::recursive_mutex> mtx1, GlobalCoolerMutex& mtx2,
           bool acquire = true);

  void lock();
  void unlock();
};

}  // namespace hictkpy
