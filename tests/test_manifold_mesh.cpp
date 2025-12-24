#include <cassert>
#include <cstddef>
#include <ranges>
#include <vector>

#include "mesh1/manifold_mesh.hpp"

namespace {

struct Empty {
  [[nodiscard]] constexpr bool operator==(const Empty&) const = default;
};

template <class Range>
std::size_t count_range(Range&& range) {
  std::size_t n = 0;
  for (auto&& _ : range) {
    (void)_;
    ++n;
  }
  return n;
}

}  // namespace

void test_manifold_mesh_single_triangle_boundary_loop() {
  using Mesh = mesh1::ManifoldMesh<Empty, Empty, Empty, Empty>;

  const Mesh mesh = Mesh::new_in(std::vector<std::vector<std::size_t>>{{0, 1, 2}});

  assert(mesh.n_vertices() == 3);
  assert(mesh.n_faces() == 1);
  assert(mesh.n_edges() == 3);
  assert(mesh.n_halfedges() == 6);

  assert(count_range(mesh.face(mesh1::FaceId{0}).halfedges()) == 3);

  std::vector<mesh1::HalfedgeId> boundary_halfedges;
  for (const auto he : mesh.halfedges()) {
    if (!mesh.he_face(he.id).valid()) {
      boundary_halfedges.push_back(he.id);
    }
  }
  assert(boundary_halfedges.size() == 3);

  {
    const auto start = boundary_halfedges[0];
    auto curr = start;
    std::size_t count = 0;
    while (true) {
      ++count;
      curr = mesh.he_next(curr);
      if (curr == start) {
        break;
      }
      assert(count < 10);
    }
    assert(count == 3);
  }

  for (const mesh1::VertexId vid : {mesh1::VertexId{0}, mesh1::VertexId{1}, mesh1::VertexId{2}}) {
    assert(count_range(mesh.vertex(vid).incoming_halfedges()) == 2);
    assert(count_range(mesh.vertex(vid).outgoing_halfedges()) == 2);
  }

  assert(count_range(mesh.edge(mesh1::EdgeId{0}).halfedges()) == 2);
}

void test_manifold_mesh_tetrahedron_closed() {
  using Mesh = mesh1::ManifoldMesh<Empty, Empty, Empty, Empty>;

  const Mesh mesh = Mesh::new_in(std::vector<std::vector<std::size_t>>{
      {0, 1, 2},
      {0, 2, 3},
      {0, 3, 1},
      {1, 3, 2},
  });

  assert(mesh.n_vertices() == 4);
  assert(mesh.n_faces() == 4);
  assert(mesh.n_edges() == 6);
  assert(mesh.n_halfedges() == 12);

  assert(count_range(mesh.halfedges() | std::views::filter([&](const auto he) {
           return !mesh.he_face(he.id).valid();
         })) == 0);

  for (std::size_t eid = 0; eid < mesh.n_edges_capacity(); ++eid) {
    assert(count_range(mesh.edge(mesh1::EdgeId{eid}).halfedges()) == 2);
  }
}
