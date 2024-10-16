// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "hictkpy/nanobind.hpp"

namespace hictkpy {

template <typename T>
struct Dynamic1DA {
 private:
  using BufferT = nanobind::ndarray<nanobind::numpy, nanobind::shape<-1>, T>;
  using VectorT = decltype(std::declval<BufferT>().view());
  nanobind::object _dtype{};
  nanobind::object _np_array{};
  BufferT _buff{};
  VectorT _vector{};

  std::int64_t _size{};
  std::int64_t _capacity{};

 public:
  explicit Dynamic1DA(std::size_t size_ = 1000);
  void push_back(T x);
  void emplace_back(T &&x);
  void resize(std::int64_t new_size);
  void grow();
  void shrink_to_fit();
  [[nodiscard]] auto operator()() -> BufferT;

 private:
  [[nodiscard]] static nanobind::object np();
};

}  // namespace hictkpy

#include "./impl/dynamic_1d_array_impl.hpp"
