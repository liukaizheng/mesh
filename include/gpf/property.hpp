#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <cmath>

#include "gpf/handles.hpp"


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

template<std::size_t N, class Mesh>
void update_edge_length_squared(gpf::EdgeHandle<Mesh, false> edge) {
    const auto [va, vb] = edge.vertices();
    const auto pa = position_span<N>(va.prop().pt);
    const auto pb = position_span<N>(vb.prop().pt);
    edge.prop().square_len = squared_distance(pa, pb);
}

template <std::size_t N, class Mesh>
void update_edge_lengths_squared(Mesh& mesh) {
  for (auto e : mesh.edges()) {
    update_edge_length_squared<N>(e);
  }
}

template<std::size_t N, class Mesh>
void update_edge_length(gpf::EdgeHandle<Mesh, false> edge) {
    const auto [va, vb] = edge.vertices();
    const auto pa = position_span<N>(va.prop().pt);
    const auto pb = position_span<N>(vb.prop().pt);
    edge.prop().len = std::sqrt(squared_distance(pa, pb));
}

template <std::size_t N, class Mesh>
void update_edge_lengths(Mesh& mesh) {
  for (auto e : mesh.edges()) {
    update_edge_length<N>(e);
  }
}

template<class Mesh>
void update_corner_angle(
  gpf::HalfedgeHandle<Mesh, false> hab,
  gpf::HalfedgeHandle<Mesh, false> hbc,
  gpf::HalfedgeHandle<Mesh, false> hca
) {
  auto lab = hab.edge().prop().len;
  auto lbc = hbc.edge().prop().len;
  auto lca = hca.edge().prop().len;
  auto q = (lab * lab + lbc * lbc - lca * lca) / (2.0 * lab * lbc);
  hca.prop().angle = std::acos(std::max(-1.0, std::min(q, 1.0)));
}

template <class Mesh>
void update_corner_angles(Mesh& mesh) {
  for (auto face : mesh.faces()) {
    auto ha = face.halfedge();
    auto hb = ha.next();
    auto hc = hb.next();
    assert(hc.next().id == ha.id);
    update_corner_angle(ha, hb, hc);
    update_corner_angle(hb, hc, ha);
    update_corner_angle(hc, ha, hb);
  }
}
template <class Mesh>
void update_vertex_angle_sum(gpf::VertexHandle<Mesh, false>& vertex) {
  double sum = 0.0;
  for (auto he : vertex.incoming_halfedges()) {
    sum += he.prop().angle;
  }
  vertex.prop().angle_sum = sum;
}

template <class Mesh>
void update_vertex_angle_sums(Mesh& mesh) {
  for (auto vertex : mesh.vertices()) {
    update_vertex_angle_sum(vertex);
  }
}
}  // namespace gpf
