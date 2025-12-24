#pragma once

#include <compare>
#include <cstddef>
#include <functional>
#include <limits>

namespace mesh1 {

inline constexpr std::size_t kInvalidIndex = std::numeric_limits<std::size_t>::max();

struct VertexId {
  std::size_t idx{kInvalidIndex};
  constexpr VertexId() = default;
  constexpr explicit VertexId(std::size_t index) : idx(index) {}
  [[nodiscard]] constexpr bool valid() const { return idx != kInvalidIndex; }
  [[nodiscard]] constexpr std::size_t index() const { return idx; }
  constexpr auto operator<=>(const VertexId&) const = default;
};

struct HalfedgeId {
  std::size_t idx{kInvalidIndex};
  constexpr HalfedgeId() = default;
  constexpr explicit HalfedgeId(std::size_t index) : idx(index) {}
  [[nodiscard]] constexpr bool valid() const { return idx != kInvalidIndex; }
  [[nodiscard]] constexpr std::size_t index() const { return idx; }
  constexpr auto operator<=>(const HalfedgeId&) const = default;
};

struct FaceId {
  std::size_t idx{kInvalidIndex};
  constexpr FaceId() = default;
  constexpr explicit FaceId(std::size_t index) : idx(index) {}
  [[nodiscard]] constexpr bool valid() const { return idx != kInvalidIndex; }
  [[nodiscard]] constexpr std::size_t index() const { return idx; }
  constexpr auto operator<=>(const FaceId&) const = default;
};

struct EdgeId {
  std::size_t idx{kInvalidIndex};
  constexpr EdgeId() = default;
  constexpr explicit EdgeId(std::size_t index) : idx(index) {}
  [[nodiscard]] constexpr bool valid() const { return idx != kInvalidIndex; }
  [[nodiscard]] constexpr std::size_t index() const { return idx; }
  constexpr auto operator<=>(const EdgeId&) const = default;
};

[[nodiscard]] constexpr HalfedgeId operator+(HalfedgeId id, std::size_t offset) {
  return HalfedgeId{id.idx + offset};
}

[[nodiscard]] constexpr HalfedgeId operator^(HalfedgeId id, std::size_t mask) {
  return HalfedgeId{id.idx ^ mask};
}

}  // namespace mesh1

namespace std {

template <>
struct hash<mesh1::VertexId> {
  std::size_t operator()(const mesh1::VertexId id) const noexcept { return id.idx; }
};

template <>
struct hash<mesh1::HalfedgeId> {
  std::size_t operator()(const mesh1::HalfedgeId id) const noexcept { return id.idx; }
};

template <>
struct hash<mesh1::FaceId> {
  std::size_t operator()(const mesh1::FaceId id) const noexcept { return id.idx; }
};

template <>
struct hash<mesh1::EdgeId> {
  std::size_t operator()(const mesh1::EdgeId id) const noexcept { return id.idx; }
};

}  // namespace std
