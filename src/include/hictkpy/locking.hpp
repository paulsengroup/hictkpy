// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef __SANITIZE_THREAD__
#define HICTKPY_USE_TSAN
#elif defined(__has_feature) && __has_feature(thread_sanitizer)
#define HICTKPY_USE_TSAN
#endif

#ifdef HICTKPY_USE_TSAN
#include <sanitizer/tsan_interface.h>
#endif

#include <Python.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <mutex>
#include <thread>

namespace hictkpy {

[[nodiscard]] constexpr bool tsan_enabled() noexcept {
#ifdef HICTKPY_USE_TSAN
  return true;
#else
  return false;
#endif
}

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

template <bool NO_LOG>
class GilScopedAcquire {
  PyGILState_STATE _state;
#ifdef HICTKPY_USE_TSAN
  static inline char _tsan_gil_proxy{};
#endif

  struct GilMutexProxyRAII {  // NOLINT(*-special-member-functions)
    explicit GilMutexProxyRAII() noexcept {
#ifdef HICTKPY_USE_TSAN
      __tsan_mutex_create(static_cast<void*>(&_tsan_gil_proxy), __tsan_mutex_write_reentrant);
#endif
    }
    ~GilMutexProxyRAII() noexcept {
#ifdef HICTKPY_USE_TSAN
      __tsan_mutex_destroy(static_cast<void*>(&_tsan_gil_proxy), __tsan_mutex_linker_init);
#endif
    }
  };

  [[nodiscard]] static PyGILState_STATE state_ensure() noexcept {
    if constexpr (!NO_LOG) {
      SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: acquiring..."), std::this_thread::get_id());
    }
    auto state = PyGILState_Ensure();
    tsan_acquire();
    if constexpr (!NO_LOG) {
      SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: acquired!"), std::this_thread::get_id());
    }
    return state;
  }

  explicit GilScopedAcquire() noexcept : _state(state_ensure()) {}

 public:
  [[nodiscard]] static auto create() { return GilScopedAcquire{}; }

  GilScopedAcquire(const GilScopedAcquire&) = delete;
  GilScopedAcquire(GilScopedAcquire&&) noexcept = default;
  ~GilScopedAcquire() noexcept {
    if constexpr (!NO_LOG) {
      SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: releasing..."), std::this_thread::get_id());
    }
    tsan_release();
    PyGILState_Release(_state);
    if constexpr (!NO_LOG) {
      SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: released!"), std::this_thread::get_id());
    }
  }

  GilScopedAcquire& operator=(const GilScopedAcquire&) = delete;
  GilScopedAcquire& operator=(GilScopedAcquire&&) noexcept = default;

  [[nodiscard]] static const GilMutexProxyRAII& try_register_with_tsan() {
    static const GilMutexProxyRAII proxy{};
    return proxy;
  }

 private:
  static void tsan_acquire() noexcept {
#ifdef HICTKPY_USE_TSAN
    __tsan_acquire(static_cast<void*>(&_tsan_gil_proxy));
#endif
  }

  static void tsan_release() noexcept {
#ifdef HICTKPY_USE_TSAN
    __tsan_release(static_cast<void*>(&_tsan_gil_proxy));
#endif
  }
};

}  // namespace hictkpy

#define HICTKPY_LOCK_COOLER_MTX_SCOPED \
  [[maybe_unused]] const auto cooler_lock = hictkpy::CoolerGlobalLock::lock();

#define HICTKPY_GIL_SCOPED_ACQUIRE \
  [[maybe_unused]] const auto hictkpy_gil_acquire_ = hictkpy::GilScopedAcquire<false>::create();

/*
#define HICTKPY_GIL_SCOPED_RELEASE \
  [[maybe_unused]] const hictkpy::GilScopedRelease hictkpy_gil_release_{};
*/
