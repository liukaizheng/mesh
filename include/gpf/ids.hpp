#pragma once

#include <compare>
#include <cstddef>
#include <functional>
#include <limits>

namespace gpf {

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

}  // namespace gpf

namespace std {

template <>
struct hash<gpf::VertexId> {
  std::size_t operator()(const gpf::VertexId id) const noexcept { return id.idx; }
};

template <>
struct hash<gpf::HalfedgeId> {
  std::size_t operator()(const gpf::HalfedgeId id) const noexcept { return id.idx; }
};

template <>
struct hash<gpf::FaceId> {
  std::size_t operator()(const gpf::FaceId id) const noexcept { return id.idx; }
};

template <>
struct hash<gpf::EdgeId> {
  std::size_t operator()(const gpf::EdgeId id) const noexcept { return id.idx; }
};

}  // namespace std
