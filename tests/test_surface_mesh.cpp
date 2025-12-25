#include <cassert>
#include <cstddef>
#include <unordered_set>
#include <vector>
#include <ranges>

#include "gpf/ids.hpp"
#include "gpf/surface_mesh.hpp"

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

void test_surface_mesh_basic() {
  using Mesh = gpf::SurfaceMesh<Empty, Empty, Empty, Empty>;

  Mesh mesh = Mesh::new_in(std::vector<std::vector<std::size_t>>{
      {0, 1, 2},
      {0, 2, 3},
      {0, 3, 1},
      {0, 4, 5},
      {0, 5, 6},
      {0, 6, 4},
  });

  assert(mesh.n_vertices() == 7);

  const auto incoming_0 = count_range(mesh.vertex(gpf::VertexId{0}).incoming_halfedges());
  const auto incomming_halfedges = mesh.vertex(gpf::VertexId{0}).incoming_halfedges() |
                                   std::views::transform([](const auto& h) { return h.id; }) |
                                   std::ranges::to<std::vector>();
  assert(incoming_0 == 6);

  const auto outgoing_0 = count_range(mesh.vertex(gpf::VertexId{0}).outgoing_halfedges());
  assert(outgoing_0 == 6);

  const auto edges_0 = count_range(mesh.vertex(gpf::VertexId{0}).edges());
  assert(edges_0 == 6);

  const auto verts_0 = count_range(mesh.vertex(gpf::VertexId{0}).vertices());
  assert(verts_0 == 6);

  const auto incoming_6 = count_range(mesh.vertex(gpf::VertexId{6}).incoming_halfedges());
  assert(incoming_6 == 2);

  const auto outgoing_6 = count_range(mesh.vertex(gpf::VertexId{6}).outgoing_halfedges());
  assert(outgoing_6 == 2);

  const auto edges_6 = count_range(mesh.vertex(gpf::VertexId{6}).edges());
  assert(edges_6 == 3);

  const auto verts_6 = count_range(mesh.vertex(gpf::VertexId{6}).vertices());
  assert(verts_6 == 3);

  const auto edge_halfedges = count_range(mesh.edge(gpf::EdgeId{0}).halfedges());
  assert(edge_halfedges == 2);

  const auto face_halfedges = count_range(mesh.face(gpf::FaceId{0}).halfedges());
  assert(face_halfedges == 3);

  const auto face_halfedges_rev = count_range(mesh.face(gpf::FaceId{0}).halfedges_reverse());
  assert(face_halfedges_rev == 3);
}

void test_surface_mesh_split_edge_and_split_face() {
  using Mesh = gpf::SurfaceMesh<Empty, Empty, Empty, Empty>;

  Mesh mesh = Mesh::new_in(std::vector<std::vector<std::size_t>>{
      {0, 1, 2},
      {0, 1, 3},
      {1, 0, 4},
      {0, 1, 5},
      {1, 0, 6},
  });

  const auto eid = mesh.e_from_vertices(gpf::VertexId{0}, gpf::VertexId{1});
  assert(eid.valid());

  const std::size_t edge_n_halfedges = count_range(mesh.edge(eid).halfedges());
  assert(edge_n_halfedges > 2);

  const auto new_vid = mesh.split_edge(eid);
  assert(new_vid.valid());

  for (auto f : mesh.faces()) {
    std::vector<gpf::HalfedgeId> hs;
    for (auto he : f.halfedges()) {
      hs.push_back(he.id);
    }
    assert(hs.size() == 4);

    for (std::size_t i = 0; i < hs.size(); ++i) {
      const auto h1 = hs[i];
      const auto h2 = hs[(i + 1) % hs.size()];
      assert(mesh.halfedge_data(h1).next == h2);
      assert(mesh.halfedge_data(h2).prev == h1);
    }

    std::unordered_set<gpf::VertexId> vset;
    for (const auto hid : hs) {
      const auto& he = mesh.halfedge_data(hid);
      assert(he.face == f.id);
      vset.insert(he.vertex);
    }
    assert(vset.size() == 4);
  }

  {
    std::vector<gpf::HalfedgeId> halfedges;
    for (auto he : mesh.edge(eid).halfedges()) {
      halfedges.push_back(he.id);
    }
    assert(halfedges.size() == edge_n_halfedges);
    for (const auto hid : halfedges) {
      assert(mesh.halfedge_data(hid).edge == eid);
    }
  }

  const auto new_eid = mesh.vertex(new_vid).halfedge().prev().edge().id;
  assert(new_eid.valid());
  {
    std::vector<gpf::HalfedgeId> halfedges;
    for (auto he : mesh.edge(new_eid).halfedges()) {
      halfedges.push_back(he.id);
    }
    assert(halfedges.size() == edge_n_halfedges);
    for (const auto hid : halfedges) {
      assert(mesh.halfedge_data(hid).edge == new_eid);
    }
  }

  {
    const gpf::FaceId fid{0};
    const auto new_hid = mesh.split_face(fid, new_vid, gpf::VertexId{2});
    const auto new_fid = mesh.he_face(new_hid);
    assert(new_fid.valid());

    std::vector<gpf::HalfedgeId> old_halfedges;
    for (auto he : mesh.face(fid).halfedges()) {
      old_halfedges.push_back(he.id);
    }
    std::vector<gpf::HalfedgeId> new_halfedges;
    for (auto he : mesh.face(new_fid).halfedges()) {
      new_halfedges.push_back(he.id);
    }

    assert(old_halfedges.size() == 3);
    assert(new_halfedges.size() == 3);

    for (std::size_t i = 0; i < old_halfedges.size(); ++i) {
      const auto h1 = old_halfedges[i];
      const auto h2 = old_halfedges[(i + 1) % old_halfedges.size()];
      assert(mesh.halfedge_data(h1).next == h2);
      assert(mesh.halfedge_data(h2).prev == h1);
      assert(mesh.halfedge_data(h1).face == fid);
    }

    for (std::size_t i = 0; i < new_halfedges.size(); ++i) {
      const auto h1 = new_halfedges[i];
      const auto h2 = new_halfedges[(i + 1) % new_halfedges.size()];
      assert(mesh.halfedge_data(h1).next == h2);
      assert(mesh.halfedge_data(h2).prev == h1);
      assert(mesh.halfedge_data(h1).face == new_fid);
    }
  }
}
