#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>

#include "gpf/handles.hpp"

namespace gpf {

template <std::size_t N, typename Mesh>
[[nodiscard]] std::span<const double, N> position_span(const VertexHandle<Mesh, false>& v) {
  return std::span<const double, N>{v.prop().pt};
}

template <std::size_t N>
[[nodiscard]] double squared_distance(std::span<const double, N> a,
                                      std::span<const double, N> b) {
  double sum = 0.0;
  for (std::size_t i = 0; i < N; ++i) {
    const double d = a[i] - b[i];
    sum += d * d;
  }
  return sum;
}

template <std::size_t N, class F, class Arg>
concept ReturnsPositionSpan = requires(F f, Arg arg) {
  { f(arg) } -> std::same_as<std::span<const double, N>>;
};

template <std::size_t N, class Mesh, class PositionSpan>
  requires ReturnsPositionSpan<N, PositionSpan, const VertexHandle<Mesh, false>&>
void update_edge_length_squared(
    gpf::EdgeHandle<Mesh, false> edge,
    PositionSpan pos_span) {
  const auto [va, vb] = edge.vertices();
  const auto pa = pos_span(va);
  const auto pb = pos_span(vb);
  edge.prop().square_len = squared_distance(pa, pb);
}

template <std::size_t N, class Mesh>
void update_edge_length_squared(gpf::EdgeHandle<Mesh, false> edge) {
  update_edge_length_squared<N>(edge, position_span<N, Mesh>);
}

template <std::size_t N, class Mesh, class PositionSpan>
  requires ReturnsPositionSpan<N, PositionSpan, const VertexHandle<Mesh, false>&>
void update_edge_lengths_squared(Mesh &mesh, PositionSpan pos_span) {
  for (auto e : mesh.edges()) {
    update_edge_length_squared<N>(e, pos_span);
  }
}

template <std::size_t N, class Mesh>
void update_edge_lengths_squared(Mesh &mesh) {
  update_edge_lengths_squared<N>(mesh, position_span<N, Mesh>);
}

template <std::size_t N, class Mesh, class PositionSpan>
  requires ReturnsPositionSpan<N, PositionSpan, const VertexHandle<Mesh, false>&>
void update_edge_length(gpf::EdgeHandle<Mesh, false> edge,
                        PositionSpan pos_span) {
  const auto [va, vb] = edge.vertices();
  const auto pa = pos_span(va);
  const auto pb = pos_span(vb);
  edge.prop().len = std::sqrt(squared_distance(pa, pb));
}

template <std::size_t N, class Mesh>
void update_edge_length(gpf::EdgeHandle<Mesh, false> edge) {
  update_edge_length<N>(edge, position_span<N, Mesh>);
}

template <std::size_t N, class Mesh, class PositionSpan>
  requires ReturnsPositionSpan<N, PositionSpan, const VertexHandle<Mesh, false>&>
void update_edge_lengths(Mesh &mesh, PositionSpan pos_span) {
  for (auto e : mesh.edges()) {
    update_edge_length<N>(e, pos_span);
  }
}

template <std::size_t N, class Mesh>
void update_edge_lengths(Mesh &mesh) {
  update_edge_lengths<N>(mesh, position_span<N, Mesh>);
}

template <class Mesh>
void update_corner_angle(gpf::HalfedgeHandle<Mesh, false> hab,
                         gpf::HalfedgeHandle<Mesh, false> hbc,
                         gpf::HalfedgeHandle<Mesh, false> hca) {
  auto lab = hab.edge().prop().len;
  auto lbc = hbc.edge().prop().len;
  auto lca = hca.edge().prop().len;
  auto q = (lab * lab + lbc * lbc - lca * lca) / (2.0 * lab * lbc);
  hca.prop().angle = std::acos(std::max(-1.0, std::min(q, 1.0)));
}

template <class Mesh> void update_corner_angles(Mesh &mesh) {
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
void update_vertex_angle_sum(gpf::VertexHandle<Mesh, false> &vertex) {
  double sum = 0.0;
  for (auto he : vertex.incoming_halfedges()) {
    sum += he.prop().angle;
  }
  vertex.prop().angle_sum = sum;
}

template <class Mesh> void update_vertex_angle_sums(Mesh &mesh) {
  for (auto vertex : mesh.vertices()) {
    update_vertex_angle_sum(vertex);
  }
}
} // namespace gpf
