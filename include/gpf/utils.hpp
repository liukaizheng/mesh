#pragma once

#include <cstddef>
#include <utility>

namespace gpf {

[[nodiscard]] constexpr std::size_t oriented_index(std::size_t idx, bool reversed) {
  return (idx << 1) | (reversed ? std::size_t{1} : std::size_t{0});
}

[[nodiscard]] constexpr std::size_t strip_orientation(std::size_t idx) { return idx >> 1; }

[[nodiscard]] constexpr std::size_t twin_index(std::size_t idx) { return idx ^ 1; }

[[nodiscard]] constexpr bool is_positive(std::size_t idx) { return (idx & 1) == 0; }

[[nodiscard]] constexpr bool is_negative(std::size_t idx) { return (idx & 1) == 1; }

[[nodiscard]] constexpr std::pair<std::size_t, bool> decode_index(std::size_t idx) {
  return {strip_orientation(idx), is_negative(idx)};
}

}  // namespace gpf
