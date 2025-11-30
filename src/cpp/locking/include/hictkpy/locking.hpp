// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// check if the compiler is invoked with -fsanitize=thread
// the logic is a bit convoluted to workaround preprocessor bugs
// when compiling with older versions of GCC
#if defined(__GNUC__) && !defined(__clang__)
#ifdef __SANITIZE_THREAD__
#define HICTKPY_USE_TSAN
#endif
#else
#if defined(__has_feature) && __has_feature(thread_sanitizer)
#define HICTKPY_USE_TSAN
#endif
#endif

#ifdef HICTKPY_USE_TSAN
#include <sanitizer/tsan_interface.h>
#endif

#include <Python.h>

#include <mutex>

#define HICTKPY_LOCK_COOLER_MTX_SCOPED \
  [[maybe_unused]] const auto cooler_lock = hictkpy::CoolerGlobalLock::lock();

#define HICTKPY_GIL_SCOPED_ACQUIRE \
  [[maybe_unused]] const auto hictkpy_gil_acquire = hictkpy::GilScopedAcquire{};

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
    UniqueLock(Mutex& mtx);  // NOLINT(*-explicit-conversions)

    UniqueLock(const UniqueLock&) = delete;
    UniqueLock(UniqueLock&&) noexcept = default;
    ~UniqueLock() noexcept;

    UniqueLock& operator=(const UniqueLock&) = delete;
    UniqueLock& operator=(UniqueLock&&) noexcept = default;

    [[nodiscard]] auto get() noexcept -> std::unique_lock<Mutex>&;
    [[nodiscard]] auto get() const noexcept -> const std::unique_lock<Mutex>&;
  };

  [[nodiscard]] static auto mtx() -> UniqueLock::Mutex&;
  [[nodiscard]] static auto lock() -> UniqueLock;
};

class GilScopedAcquire {
  PyGILState_STATE _state;
#ifdef HICTKPY_USE_TSAN
  static inline char _tsan_gil_proxy{};
#endif
  bool _no_log;

  struct GilMutexProxyRAII {
    explicit GilMutexProxyRAII() noexcept;
    GilMutexProxyRAII(const GilMutexProxyRAII&) = delete;
    GilMutexProxyRAII(GilMutexProxyRAII&&) noexcept = default;
    ~GilMutexProxyRAII() noexcept;
    GilMutexProxyRAII& operator=(const GilMutexProxyRAII&) = delete;
    GilMutexProxyRAII& operator=(GilMutexProxyRAII&&) noexcept = default;
  };

  [[nodiscard]] static PyGILState_STATE state_ensure(bool no_log) noexcept;

 public:
  explicit GilScopedAcquire(bool no_log = false) noexcept;
  GilScopedAcquire(const GilScopedAcquire&) = delete;
  GilScopedAcquire(GilScopedAcquire&&) noexcept = default;
  ~GilScopedAcquire() noexcept;
  GilScopedAcquire& operator=(const GilScopedAcquire&) = delete;
  GilScopedAcquire& operator=(GilScopedAcquire&&) noexcept = default;

  [[nodiscard]] static auto try_register_with_tsan() -> const GilMutexProxyRAII&;

 private:
  static void tsan_acquire() noexcept;
  static void tsan_release() noexcept;
};

}  // namespace hictkpy
