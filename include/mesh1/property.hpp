#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <type_traits>

#include "mesh1/detail.hpp"

namespace mesh1 {

template <std::size_t N>
[[nodiscard]] std::span<const double, N> position_span(const std::array<double, N>& p) {
  return std::span<const double, N>{p};
}

template <std::size_t N, class P>
requires requires(const P& p) {
  { p.position() } -> std::same_as<const std::array<double, N>&>;
}
[[nodiscard]] std::span<const double, N> position_span(const P& p) {
  return position_span<N>(p.position());
}

template <std::size_t N, class VertexData>
requires requires(const VertexData& v) { v.property; }
[[nodiscard]] std::span<const double, N> position_span(const VertexData& v) {
  return position_span<N>(v.property);
}

template <std::size_t N>
[[nodiscard]] double squared_distance(std::span<const double, N> a, std::span<const double, N> b) {
  double sum = 0.0;
  for (std::size_t i = 0; i < N; ++i) {
    const double d = a[i] - b[i];
    sum += d * d;
  }
  return sum;
}

template <std::size_t N, class Mesh>
void update_edge_lengths_squared_in_edge_data(Mesh& mesh) {
  for (auto e : mesh.edges()) {
    const auto [va, vb] = e.vertices();
    if (!va.id.valid() || !vb.id.valid()) {
      continue;
    }

    const auto pa = position_span<N>(va.data());
    const auto pb = position_span<N>(vb.data());
    auto& prop = detail::property_ref(e.data());
    detail::edge_length_squared_mut(prop) = squared_distance<N>(pa, pb);
  }
}

template <std::size_t N, class Mesh>
void update_edge_lengths_in_edge_data(Mesh& mesh) {
  for (auto e : mesh.edges()) {
    const auto [va, vb] = e.vertices();
    if (!va.id.valid() || !vb.id.valid()) {
      continue;
    }

    const auto pa = position_span<N>(va.data());
    const auto pb = position_span<N>(vb.data());
    auto& prop = detail::property_ref(e.data());
    detail::edge_length_mut(prop) = std::sqrt(squared_distance<N>(pa, pb));
  }
}

}  // namespace mesh1

