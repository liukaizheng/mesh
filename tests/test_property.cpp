#include <cassert>
#include <cmath>
#include <array>
#include <cstddef>
#include <vector>

#include "gpf/property.hpp"
#include "gpf/surface_mesh.hpp"

namespace {

struct Empty {
  [[nodiscard]] constexpr bool operator==(const Empty&) const = default;
};

struct EdgeProp {
  double square_len = 0.0;
  double len = 0.0;
};

}  // namespace

void test_property_edge_length_updates() {
  using Mesh = gpf::SurfaceMesh<std::array<double, 3>, Empty, EdgeProp, Empty>;

  Mesh mesh = Mesh::new_in(std::vector<std::vector<std::size_t>>{{0, 1, 2}});

  for (auto v : mesh.vertices()) {
    const auto id = v.id.idx;
    if (id == 0) {
      v.data().property = {0.0, 0.0, 0.0};
    } else if (id == 1) {
      v.data().property = {1.0, 0.0, 0.0};
    } else if (id == 2) {
      v.data().property = {0.0, 2.0, 0.0};
    } else {
      assert(false);
    }
  }

  gpf::update_edge_lengths_squared_in_edge_data<3>(mesh);
  gpf::update_edge_lengths_in_edge_data<3>(mesh);

  for (auto e : mesh.edges()) {
    const auto [va, vb] = e.vertices();
    const auto pa = va.data().property;
    const auto pb = vb.data().property;

    const double dx = pa[0] - pb[0];
    const double dy = pa[1] - pb[1];
    const double dz = pa[2] - pb[2];
    const double expected = dx * dx + dy * dy + dz * dz;
    const double expected_len = std::sqrt(expected);

    assert(std::abs(e.data().property.square_len - expected) < 1e-12);
    assert(std::abs(e.data().property.len - expected_len) < 1e-12);
  }
}
