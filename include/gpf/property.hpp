#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <span>

#include "gpf/detail.hpp"

namespace gpf {

template <std::size_t N>
[[nodiscard]] std::span<const double, N> position_span(const std::array<double, N>& p) {
  return std::span<const double, N>{p};
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
void update_edge_lengths_squared(Mesh& mesh) {
  for (auto e : mesh.edges()) {
    const auto [va, vb] = e.vertices();
    const auto pa = position_span<N>(va.prop().pt);
    const auto pb = position_span<N>(vb.prop().pt);
    e.prop().square_len = squared_distance(pa, pb);
  }
}

template <std::size_t N, class Mesh>
void update_edge_lengths(Mesh& mesh) {
  for (auto e : mesh.edges()) {
    const auto [va, vb] = e.vertices();
    const auto pa = position_span<N>(va.prop().pt);
    const auto pb = position_span<N>(vb.prop().pt);
    e.prop().len = std::sqrt(squared_distance<N>(pa, pb));
  }
}
}  // namespace gpf
