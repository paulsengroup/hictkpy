// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "hictkpy/common.hpp"
#include "hictkpy/nanobind.hpp"

namespace hictkpy {

template <typename T>
inline Dynamic1DA<T>::Dynamic1DA(std::size_t size_)
    : _dtype(np().attr("dtype")(map_type_to_dtype<T>())),
      _np_array(
          np().attr("empty")(static_cast<std::int64_t>(size_), nanobind::arg("dtype") = _dtype)),
      _buff(nanobind::cast<BufferT>(_np_array)),
      _vector(_buff.view()),
      _capacity(static_cast<std::int64_t>(size_)) {}

template <typename T>
inline void Dynamic1DA<T>::push_back(T x) {
  if (_capacity == _size) {
    grow();
  }
  _vector(_size++) = x;
}

template <typename T>
inline void Dynamic1DA<T>::emplace_back(T &&x) {
  if (_capacity == _size) {
    grow();
  }
  _vector(_size++) = std::move(x);
}

template <typename T>
inline void Dynamic1DA<T>::resize(std::int64_t new_size) {
  if (_capacity == new_size) {
    return;
  }
  auto new_array = np().attr("empty")(new_size, nanobind::arg("dtype") = _dtype);
  auto new_buff = nanobind::cast<BufferT>(new_array);
  auto new_vector = new_buff.view();

  _capacity = new_size;
  _size = std::min(_capacity, _size);
  std::copy(_vector.data(), _vector.data() + static_cast<std::size_t>(_size), new_vector.data());

  std::swap(new_array, _np_array);
  std::swap(new_buff, _buff);
  std::swap(new_vector, _vector);
}

template <typename T>
inline void Dynamic1DA<T>::grow() {
  resize(static_cast<std::int64_t>(_buff.size() * 2));
}

template <typename T>
inline void Dynamic1DA<T>::shrink_to_fit() {
  resize(_size);
}

template <typename T>
[[nodiscard]] auto Dynamic1DA<T>::operator()() -> BufferT {
  shrink_to_fit();
  return _buff;
}

template <typename T>
[[nodiscard]] nanobind::object Dynamic1DA<T>::np() {
  return nanobind::module_::import_("numpy");
}

}  // namespace hictkpy
