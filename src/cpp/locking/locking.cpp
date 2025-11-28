// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictkpy/locking.hpp"

#include <Python.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <mutex>
#include <thread>

namespace hictkpy {
CoolerGlobalLock::UniqueLock::UniqueLock(Mutex& mtx) : _lck(mtx) {
  SPDLOG_TRACE(FMT_STRING("[tid={}]: CoolerGlobalLock({}): locked!"), std::this_thread::get_id(),
               fmt::ptr(&mtx));
}

CoolerGlobalLock::UniqueLock::~UniqueLock() noexcept {
  if (const auto* mtx = _lck.mutex(); !!mtx) {
    SPDLOG_TRACE(FMT_STRING("[tid={}]: CoolerGlobalLock({}): unlocking..."),
                 std::this_thread::get_id(), fmt::ptr(mtx));
    _lck.unlock();
    SPDLOG_TRACE(FMT_STRING("[tid={}]: CoolerGlobalLock({}): unlocked!"),
                 std::this_thread::get_id(), fmt::ptr(mtx));
  }
}

auto CoolerGlobalLock::UniqueLock::get() noexcept -> std::unique_lock<Mutex>& { return _lck; }
auto CoolerGlobalLock::UniqueLock::get() const noexcept -> const std::unique_lock<Mutex>& {
  return _lck;
}

auto CoolerGlobalLock::mtx() -> UniqueLock::Mutex& {
  static UniqueLock::Mutex mtx_;
  return mtx_;
}

auto CoolerGlobalLock::lock() -> UniqueLock {
  SPDLOG_TRACE(FMT_STRING("[tid={}]: CoolerGlobalLock({}): locking..."), std::this_thread::get_id(),
               fmt::ptr(&mtx()));
  return {mtx()};
}

GilScopedAcquire::GilMutexProxyRAII::GilMutexProxyRAII() noexcept {
#ifdef HICTKPY_USE_TSAN
  __tsan_mutex_create(static_cast<void*>(&_tsan_gil_proxy), __tsan_mutex_write_reentrant);
#endif
}
GilScopedAcquire::GilMutexProxyRAII::~GilMutexProxyRAII() noexcept {
#ifdef HICTKPY_USE_TSAN
  __tsan_mutex_destroy(static_cast<void*>(&_tsan_gil_proxy), __tsan_mutex_linker_init);
#endif
}

PyGILState_STATE GilScopedAcquire::state_ensure(bool no_log) noexcept {
  if (!no_log) {
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: acquiring..."), std::this_thread::get_id());
  }
  auto state = PyGILState_Ensure();
  tsan_acquire();
  if (!no_log) {
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: acquired!"), std::this_thread::get_id());
  }
  return state;
}

GilScopedAcquire::GilScopedAcquire(bool no_log) noexcept
    : _state(state_ensure(no_log)), _no_log(no_log) {}

GilScopedAcquire::~GilScopedAcquire() noexcept {
  if (!_no_log) {
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: releasing..."), std::this_thread::get_id());
  }
  tsan_release();
  PyGILState_Release(_state);
  if (!_no_log) {
    SPDLOG_TRACE(FMT_STRING("[tid={}]: GIL: released!"), std::this_thread::get_id());
  }
}

auto GilScopedAcquire::try_register_with_tsan() -> const GilMutexProxyRAII& {
  static const GilMutexProxyRAII proxy{};
  return proxy;
}

void GilScopedAcquire::tsan_acquire() noexcept {
#ifdef HICTKPY_USE_TSAN
  __tsan_acquire(static_cast<void*>(&_tsan_gil_proxy));
#endif
}

void GilScopedAcquire::tsan_release() noexcept {
#ifdef HICTKPY_USE_TSAN
  __tsan_release(static_cast<void*>(&_tsan_gil_proxy));
#endif
}
}  // namespace hictkpy
