#include <algorithm>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <ranges>
#include <vector>

#include "gpf/ids.hpp"
#include "gpf/manifold_mesh.hpp"

namespace {

struct Empty {
  [[nodiscard]] constexpr bool operator==(const Empty &) const = default;
};

template <class Range> std::size_t count_range(Range &&range) {
  std::size_t n = 0;
  for (auto &&_ : range) {
    (void)_;
    ++n;
  }
  return n;
}

} // namespace

void test_manifold_mesh_single_triangle_boundary_loop() {
  using Mesh = gpf::ManifoldMesh<Empty, Empty, Empty, Empty>;

  const Mesh mesh =
      Mesh::new_in(std::vector<std::vector<std::size_t>>{{0, 1, 2}});

  assert(mesh.n_vertices() == 3);
  assert(mesh.n_faces() == 1);
  assert(mesh.n_edges() == 3);
  assert(mesh.n_halfedges() == 6);

  assert(count_range(mesh.face(gpf::FaceId{0}).halfedges()) == 3);

  std::vector<gpf::HalfedgeId> boundary_halfedges;
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

  for (const gpf::VertexId vid :
       {gpf::VertexId{0}, gpf::VertexId{1}, gpf::VertexId{2}}) {
    assert(count_range(mesh.vertex(vid).incoming_halfedges()) == 2);
    assert(count_range(mesh.vertex(vid).outgoing_halfedges()) == 2);
  }

  assert(count_range(mesh.edge(gpf::EdgeId{0}).halfedges()) == 2);
}

void test_manifold_vertex_halfedges() {
  using Mesh = gpf::ManifoldMesh<Empty, Empty, Empty, Empty>;

  const Mesh mesh = Mesh::new_in(std::vector<std::vector<std::size_t>>{
      {0, 2, 3},
      {0, 1, 2},
      {0, 3, 4}
  });

  const auto incoming_vertices = mesh.vertex(gpf::VertexId{0}).incoming_halfedges() |
                                 std::views::transform([](const auto he) { return he.from().id; }) |
                                 std::ranges::to<std::vector>();
  const auto outgoing_vertices = mesh.vertex(gpf::VertexId{0}).outgoing_halfedges() |
                                 std::views::transform([](const auto he) { return he.to().id; }) |
                                 std::ranges::to<std::vector>();
  assert(incoming_vertices.front() == gpf::VertexId{1ul});
  assert(outgoing_vertices.front() == gpf::VertexId{4ul});
}

void test_manifold_mesh_tetrahedron_closed() {
  using Mesh = gpf::ManifoldMesh<Empty, Empty, Empty, Empty>;

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

  auto vertices = mesh.vertex(gpf::VertexId{0}).outgoing_halfedges() |
                  std::views::transform([](const auto v) { return v.to().id; }) |
                  std::ranges::to<std::vector>();
  auto v_it = std::ranges::find(vertices, 2ull, &gpf::VertexId::idx);
  assert(v_it != std::ranges::end(vertices));
  std::ranges::rotate(vertices, v_it);
  assert((vertices == std::vector<gpf::VertexId>{gpf::VertexId{2ull}, gpf::VertexId{3ull}, gpf::VertexId{1ull}}));

  assert(count_range(mesh.halfedges() | std::views::filter([&](const auto he) {
                       return !mesh.he_face(he.id).valid();
                     })) == 0);

  for (std::size_t eid = 0; eid < mesh.n_edges_capacity(); ++eid) {
    assert(count_range(mesh.edge(gpf::EdgeId{eid}).halfedges()) == 2);
  }
}

void test_split_face_into_triangles() {
  using Mesh = gpf::ManifoldMesh<Empty, Empty, Empty, Empty>;

  Mesh mesh =
      Mesh::new_in(std::vector<std::vector<std::size_t>>{{0, 1, 2}});
  mesh.new_vertices(1);
  auto triangles = std::array<std::size_t, 9>{0, 1, 3, 1, 2, 3, 2, 0, 3} |
    std::views::transform([](auto v) { return gpf::VertexId{v}; }) |
    std::ranges::to<std::vector>();
  mesh.split_face_into_triangles(gpf::FaceId{0ul}, triangles);

  auto vertices = mesh.vertex(gpf::VertexId{3ul}).outgoing_halfedges() |
                  std::views::transform([](const auto he) { return he.to().id; }) |
                  std::ranges::to<std::vector>();
  assert(vertices.size() == 3);
}
