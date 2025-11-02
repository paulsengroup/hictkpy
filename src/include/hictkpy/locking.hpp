// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Python.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <mutex>
#include <thread>

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
      SPDLOG_TRACE(FMT_STRING("[tid={}]: CoolerGlobalLock({}): locked!"),
                   std::this_thread::get_id(), fmt::ptr(&mtx));
    }

    UniqueLock(const UniqueLock&) = delete;
    UniqueLock(UniqueLock&&) noexcept = default;
    ~UniqueLock() noexcept {
      if (const auto* mtx = _lck.mutex(); !!mtx) {
        SPDLOG_TRACE(FMT_STRING("[tid={}]: CoolerGlobalLock({}): unlocking..."),
                     std::this_thread::get_id(), fmt::ptr(mtx));
        _lck.unlock();
        SPDLOG_TRACE(FMT_STRING("[tid={}]: CoolerGlobalLock({}): unlocked!"),
                     std::this_thread::get_id(), fmt::ptr(mtx));
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
    SPDLOG_TRACE(FMT_STRING("[tid={}]: CoolerGlobalLock({}): locking..."),
                 std::this_thread::get_id(), fmt::ptr(&mtx()));
    return {mtx()};
  }
};

class GilScopedAcquire {
  // This mutex is likely redundant. Locking the GIL should be enough.
  // The mutex is mostly here to make TSAN happy.
  using Mutex = std::recursive_mutex;
  std::unique_lock<Mutex> _lck{};
  PyGILState_STATE _state;

  static inline Mutex _mtx{};

  [[nodiscard]] static PyGILState_STATE init() noexcept {
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: acquiring..."), std::this_thread::get_id());
    auto state = PyGILState_Ensure();
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: acquired!"), std::this_thread::get_id());
    return state;
  }

 public:
  explicit GilScopedAcquire() noexcept : _lck(_mtx), _state(init()) {}
  GilScopedAcquire(const GilScopedAcquire&) = delete;
  GilScopedAcquire(GilScopedAcquire&&) noexcept = default;
  ~GilScopedAcquire() noexcept {
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: releasing..."), std::this_thread::get_id());
    PyGILState_Release(_state);
    _lck.unlock();
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: released!"), std::this_thread::get_id());
  }

  GilScopedAcquire& operator=(const GilScopedAcquire&) = delete;
  GilScopedAcquire& operator=(GilScopedAcquire&&) noexcept = default;
};

class GilScopedRelease {
  PyThreadState* _state{};

  [[nodiscard]] static PyThreadState* init() noexcept {
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: releasing..."), std::this_thread::get_id());
    auto* state = PyEval_SaveThread();
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: released!"), std::this_thread::get_id());
    return state;
  }

 public:
  explicit GilScopedRelease() noexcept : _state(init()) {}
  GilScopedRelease(const GilScopedRelease&) = delete;
  GilScopedRelease(GilScopedRelease&&) noexcept = default;
  ~GilScopedRelease() noexcept {
    if (_state) {
      SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: acquiring..."), std::this_thread::get_id());
      PyEval_RestoreThread(_state);
      SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: acquired!"), std::this_thread::get_id());
    }
  }

  GilScopedRelease& operator=(const GilScopedRelease&) = delete;
  GilScopedRelease& operator=(GilScopedRelease&&) noexcept = default;
};

}  // namespace hictkpy

#define HICTKPY_LOCK_COOLER_MTX_SCOPED \
  [[maybe_unused]] const auto cooler_lock = hictkpy::CoolerGlobalLock::lock();

#define HICTKPY_GIL_SCOPED_ACQUIRE \
  [[maybe_unused]] const hictkpy::GilScopedAcquire hictkpy_gil_acquire_{};

#define HICTKPY_GIL_SCOPED_RELEASE \
  [[maybe_unused]] const hictkpy::GilScopedRelease hictkpy_gil_release_{};
