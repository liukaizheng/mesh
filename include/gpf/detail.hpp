#pragma once

#include <cstddef>
#include <type_traits>
#include <utility>

#include "gpf/ids.hpp"

namespace gpf::detail {

[[nodiscard]] constexpr std::pair<VertexId, VertexId> ordered_pair(VertexId a, VertexId b) {
  if (a.idx < b.idx) {
    return {a, b};
  }
  return {b, a};
}

struct PairHash {
  template <class A, class B>
  std::size_t operator()(const std::pair<A, B>& p) const noexcept {
    const std::size_t h1 = std::hash<A>{}(p.first);
    const std::size_t h2 = std::hash<B>{}(p.second);
    // https://stackoverflow.com/a/2595226
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ull + (h1 << 6) + (h1 >> 2));
  }
};

template <class T>
constexpr bool has_member_property_v = requires(T t) { t.property; };

template <class EdgeData>
decltype(auto) property_ref(EdgeData& data) {
  if constexpr (has_member_property_v<EdgeData>) {
    return (data.property);
  } else {
    return (data);
  }
}

template <class EdgeData>
decltype(auto) property_ref(const EdgeData& data) {
  if constexpr (has_member_property_v<const EdgeData>) {
    return (data.property);
  } else {
    return (data);
  }
}

template <class T>
constexpr bool has_member_square_len_v = requires(T t) { t.square_len; };

template <class T>
constexpr bool has_member_len_v = requires(T t) { t.len; };

template <class P>
double& edge_length_squared_mut(P& p) {
  if constexpr (requires { p.edge_length_squared_mut(); }) {
    return p.edge_length_squared_mut();
  } else if constexpr (has_member_square_len_v<P>) {
    return p.square_len;
  } else {
    static_assert(sizeof(P) == 0,
                  "Edge property must provide `edge_length_squared_mut()` or a `square_len` "
                  "member.");
  }
}

template <class P>
double& edge_length_mut(P& p) {
  if constexpr (requires { p.edge_length_mut(); }) {
    return p.edge_length_mut();
  } else if constexpr (has_member_len_v<P>) {
    return p.len;
  } else {
    static_assert(sizeof(P) == 0,
                  "Edge property must provide `edge_length_mut()` or a `len` member.");
  }
}

}  // namespace gpf::detail
