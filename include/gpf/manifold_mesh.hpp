#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <ranges>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gpf/detail.hpp"
#include "gpf/handles.hpp"
#include "gpf/ids.hpp"

namespace gpf {

template <class VertexProp, class HalfedgeProp, class EdgeProp, class FaceProp>
class ManifoldMesh {
 public:
  struct VertexData {
    HalfedgeId halfedge{};
    VertexProp property{};
    [[nodiscard]] bool valid() const { return halfedge.valid(); }
  };

  struct HalfedgeData {
    VertexId vertex{};
    HalfedgeId next{};
    HalfedgeId prev{};
    FaceId face{};
    HalfedgeProp property{};
    [[nodiscard]] bool valid() const { return vertex.valid(); }
  };

  struct FaceData {
    HalfedgeId halfedge{};
    FaceProp property{};
    [[nodiscard]] bool valid() const { return halfedge.valid(); }
  };

  using EdgeData = EdgeProp;

  using Self = ManifoldMesh<VertexProp, HalfedgeProp, EdgeProp, FaceProp>;
  using Vertex = VertexHandle<Self, false>;
  using Halfedge = HalfedgeHandle<Self, false>;
  using Face = FaceHandle<Self, false>;
  using Edge = EdgeHandle<Self, false>;
  using ConstVertex = VertexHandle<Self, true>;
  using ConstHalfedge = HalfedgeHandle<Self, true>;
  using ConstFace = FaceHandle<Self, true>;
  using ConstEdge = EdgeHandle<Self, true>;

  ManifoldMesh() = default;

  template <class PolygonRange>
  static Self new_in(PolygonRange&& polygons) {
    Self mesh;
    mesh.build(std::forward<PolygonRange>(polygons));
    return mesh;
  }

  [[nodiscard]] std::size_t n_vertices() const { return n_vertices_; }
  [[nodiscard]] std::size_t n_halfedges() const { return n_halfedges_; }
  [[nodiscard]] std::size_t n_faces() const { return n_faces_; }
  [[nodiscard]] std::size_t n_edges() const { return n_halfedges_ >> 1; }

  [[nodiscard]] std::size_t n_vertices_capacity() const { return vertices_.size(); }
  [[nodiscard]] std::size_t n_halfedges_capacity() const { return halfedges_.size(); }
  [[nodiscard]] std::size_t n_faces_capacity() const { return faces_.size(); }
  [[nodiscard]] std::size_t n_edges_capacity() const { return edges_.size(); }

  [[nodiscard]] const VertexData& vertex_data(VertexId vid) const { return vertices_[vid.idx]; }
  [[nodiscard]] VertexData& vertex_data(VertexId vid) { return vertices_[vid.idx]; }

  [[nodiscard]] const HalfedgeData& halfedge_data(HalfedgeId hid) const { return halfedges_[hid.idx]; }
  [[nodiscard]] HalfedgeData& halfedge_data(HalfedgeId hid) { return halfedges_[hid.idx]; }

  [[nodiscard]] const FaceData& face_data(FaceId fid) const { return faces_[fid.idx]; }
  [[nodiscard]] FaceData& face_data(FaceId fid) { return faces_[fid.idx]; }

  [[nodiscard]] const EdgeData& edge_data(EdgeId eid) const { return edges_[eid.idx]; }
  [[nodiscard]] EdgeData& edge_data(EdgeId eid) { return edges_[eid.idx]; }

  [[nodiscard]] const VertexProp& vertex_prop(VertexId vid) const { return vertex_data(vid).property; }
  [[nodiscard]] VertexProp& vertex_prop(VertexId vid) { return vertex_data(vid).property; }

  [[nodiscard]] const HalfedgeProp& halfedge_prop(HalfedgeId hid) const { return halfedge_data(hid).property; }
  [[nodiscard]] HalfedgeProp& halfedge_prop(HalfedgeId hid) { return halfedge_data(hid).property; }

  [[nodiscard]] const FaceProp& face_prop(FaceId fid) const { return face_data(fid).property; }
  [[nodiscard]] FaceProp& face_prop(FaceId fid) { return face_data(fid).property; }

  [[nodiscard]] const EdgeProp& edge_prop(EdgeId eid) const { return edge_data(eid); }
  [[nodiscard]] EdgeProp& edge_prop(EdgeId eid) { return edge_data(eid); }

  [[nodiscard]] Vertex vertex(VertexId vid) { return Vertex{vid, this}; }
  [[nodiscard]] ConstVertex vertex(VertexId vid) const { return ConstVertex{vid, this}; }
  [[nodiscard]] Halfedge halfedge(HalfedgeId hid) { return Halfedge{hid, this}; }
  [[nodiscard]] ConstHalfedge halfedge(HalfedgeId hid) const { return ConstHalfedge{hid, this}; }
  [[nodiscard]] Face face(FaceId fid) { return Face{fid, this}; }
  [[nodiscard]] ConstFace face(FaceId fid) const { return ConstFace{fid, this}; }
  [[nodiscard]] Edge edge(EdgeId eid) { return Edge{eid, this}; }
  [[nodiscard]] ConstEdge edge(EdgeId eid) const { return ConstEdge{eid, this}; }

  [[nodiscard]] auto vertices() {
    return std::views::iota(std::size_t{0}, vertices_.size()) |
           std::views::filter([this](const std::size_t i) { return vertices_[i].valid(); }) |
           std::views::transform([this](const std::size_t i) { return Vertex{VertexId{i}, this}; });
  }

  [[nodiscard]] auto vertices() const {
    return std::views::iota(std::size_t{0}, vertices_.size()) |
           std::views::filter([this](const std::size_t i) { return vertices_[i].valid(); }) |
           std::views::transform(
               [this](const std::size_t i) { return ConstVertex{VertexId{i}, this}; });
  }

  [[nodiscard]] auto halfedges() {
    return std::views::iota(std::size_t{0}, halfedges_.size()) |
           std::views::filter([this](const std::size_t i) { return halfedges_[i].valid(); }) |
           std::views::transform(
               [this](const std::size_t i) { return Halfedge{HalfedgeId{i}, this}; });
  }

  [[nodiscard]] auto halfedges() const {
    return std::views::iota(std::size_t{0}, halfedges_.size()) |
           std::views::filter([this](const std::size_t i) { return halfedges_[i].valid(); }) |
           std::views::transform(
               [this](const std::size_t i) { return ConstHalfedge{HalfedgeId{i}, this}; });
  }

  [[nodiscard]] auto faces() {
    return std::views::iota(std::size_t{0}, faces_.size()) |
           std::views::filter([this](const std::size_t i) { return faces_[i].valid(); }) |
           std::views::transform([this](const std::size_t i) { return Face{FaceId{i}, this}; });
  }

  [[nodiscard]] auto faces() const {
    return std::views::iota(std::size_t{0}, faces_.size()) |
           std::views::filter([this](const std::size_t i) { return faces_[i].valid(); }) |
           std::views::transform([this](const std::size_t i) { return ConstFace{FaceId{i}, this}; });
  }

  [[nodiscard]] auto edges() {
    return std::views::iota(std::size_t{0}, edges_.size()) |
           std::views::filter([this](const std::size_t i) { return e_is_valid(EdgeId{i}); }) |
           std::views::transform([this](const std::size_t i) { return Edge{EdgeId{i}, this}; });
  }

  [[nodiscard]] auto edges() const {
    return std::views::iota(std::size_t{0}, edges_.size()) |
           std::views::filter([this](const std::size_t i) { return e_is_valid(EdgeId{i}); }) |
           std::views::transform([this](const std::size_t i) { return ConstEdge{EdgeId{i}, this}; });
  }

  [[nodiscard]] HalfedgeId he_twin(HalfedgeId hid) const { return hid ^ 1; }

  [[nodiscard]] bool he_is_boundary(HalfedgeId hid) const { return !he_face(hid).valid(); }

  [[nodiscard]] HalfedgeId v_halfedge(VertexId vid) const { return vertex_data(vid).halfedge; }

  [[nodiscard]] VertexId he_from(HalfedgeId hid) const { return halfedges_[halfedges_[hid.idx].prev.idx].vertex; }
  [[nodiscard]] VertexId he_to(HalfedgeId hid) const { return halfedges_[hid.idx].vertex; }
  [[nodiscard]] std::array<VertexId, 2> he_vertices(HalfedgeId hid) const {
    return {he_from(hid), he_to(hid)};
  }

  [[nodiscard]] HalfedgeId he_prev(HalfedgeId hid) const { return halfedges_[hid.idx].prev; }
  [[nodiscard]] HalfedgeId he_next(HalfedgeId hid) const { return halfedges_[hid.idx].next; }
  [[nodiscard]] FaceId he_face(HalfedgeId hid) const { return halfedges_[hid.idx].face; }

  [[nodiscard]] EdgeId he_edge(HalfedgeId hid) const { return EdgeId{hid.idx >> 1}; }
  [[nodiscard]] HalfedgeId he_sibling(HalfedgeId hid) const { return he_twin(hid); }
  [[nodiscard]] HalfedgeId he_incoming_next(HalfedgeId hid) const { return he_prev(he_twin(hid)); }

  [[nodiscard]] HalfedgeId e_halfedge(EdgeId eid) const { return HalfedgeId{eid.idx << 1}; }

  [[nodiscard]] EdgeId e_from_vertices(VertexId va, VertexId vb) const {
    for (const auto e : vertex(va).edges()) {
      const auto he = e.halfedge();
      const auto v1 = he.to().id;
      if (v1 == va) {
        if (he.from().id == vb) {
          return e.id;
        }
      } else if (v1 == vb) {
        if (he.from().id == va) {
          return e.id;
        }
      }
    }
    return EdgeId{};
  }

  [[nodiscard]] std::array<VertexId, 2> e_vertices(EdgeId eid) const {
    const HalfedgeId hid = e_halfedge(eid);
    return {he_from(hid), he_to(hid)};
  }

  [[nodiscard]] HalfedgeId he_from_oppo_vertex(FaceId fid, VertexId vid) const {
    for (const auto he : face(fid).halfedges()) {
      if (he.data().vertex == vid) {
        return he.prev().id;
      }
    }
    assert(false && "vertex not found in face");
    return HalfedgeId{};
  }

  [[nodiscard]] HalfedgeId f_halfedge(FaceId fid) const { return faces_[fid.idx].halfedge; }

  void connect_halfedges(HalfedgeId hid1, HalfedgeId hid2) {
    halfedges_[hid1.idx].next = hid2;
    halfedges_[hid2.idx].prev = hid1;
  }

  VertexId new_vertices(std::size_t n) {
    const VertexId ret{vertices_.size()};
    vertices_.resize(ret.idx + n);
    n_vertices_ += n;
    return ret;
  }

  FaceId new_face() {
    const FaceId ret{faces_.size()};
    faces_.resize(ret.idx + 1);
    n_faces_ += 1;
    return ret;
  }

  HalfedgeId new_edge() {
    assert(halfedges_.size() % 2 == 0);
    const HalfedgeId hid{halfedges_.size()};
    halfedges_.resize(hid.idx + 2);
    n_halfedges_ += 2;
    assert(hid.valid());
    assert((hid.idx & 1) == 0);

    const EdgeId eid{hid.idx >> 1};
    if (edges_.size() < eid.idx + 1) {
      edges_.resize(eid.idx + 1);
    }

    return hid;
  }

  void remove_vertex(VertexId vid) {
    if (!vid.valid()) {
      return;
    }
    vertex_data(vid).halfedge = HalfedgeId{};
    n_vertices_ = (n_vertices_ == 0) ? 0 : (n_vertices_ - 1);
  }

  void remove_edge(EdgeId eid) {
    if (!eid.valid()) {
      return;
    }
    const HalfedgeId hid = e_halfedge(eid);
    const HalfedgeId twin_hid = he_twin(hid);

    halfedge_data(hid).vertex = VertexId{};
    halfedge_data(twin_hid).vertex = VertexId{};
    edges_[eid.idx] = EdgeProp{};

    n_halfedges_ = (n_halfedges_ < 2) ? 0 : (n_halfedges_ - 2);
  }

  void remove_face(FaceId fid) {
    if (!fid.valid()) {
      return;
    }

    std::vector<HalfedgeId> face_halfedges;
    for (const auto he : face(fid).halfedges()) {
      face_halfedges.push_back(he.id);
    }

    for (std::size_t i = 0; i < face_halfedges.size(); ++i) {
      const HalfedgeId curr_hid = face_halfedges[i];
      const HalfedgeId next_hid = face_halfedges[(i + 1) % face_halfedges.size()];

      halfedge_data(curr_hid).face = FaceId{};

      const HalfedgeId rev_next_hid = he_twin(curr_hid);
      const HalfedgeId rev_curr_hid = he_twin(next_hid);

      const bool rev_next_boundary = he_is_boundary(rev_next_hid);
      const bool rev_curr_boundary = he_is_boundary(rev_curr_hid);

      if (rev_next_boundary && rev_curr_boundary) {
        const VertexId vid = he_to(curr_hid);
        const HalfedgeId vh = vid.valid() ? v_halfedge(vid) : HalfedgeId{};

        if (vh == rev_next_hid && he_prev(rev_next_hid) == rev_curr_hid) {
          remove_vertex(vid);
        } else {
          const HalfedgeId prev_rev_next_hid = he_prev(rev_next_hid);
          const HalfedgeId next_rev_curr_hid = he_next(rev_curr_hid);
          connect_halfedges(prev_rev_next_hid, next_rev_curr_hid);
          if (vid.valid()) {
            vertex_data(vid).halfedge = next_rev_curr_hid;
          }
        }
      } else if (rev_next_boundary && !rev_curr_boundary) {
        connect_halfedges(he_prev(rev_next_hid), next_hid);
        const VertexId vid = he_to(curr_hid);
        if (vid.valid()) {
          vertex_data(vid).halfedge = next_hid;
        }
      } else if (!rev_next_boundary && rev_curr_boundary) {
        const HalfedgeId next_rev_curr_hid = he_next(rev_curr_hid);
        connect_halfedges(curr_hid, next_rev_curr_hid);
        const VertexId vid = he_to(curr_hid);
        if (vid.valid()) {
          vertex_data(vid).halfedge = next_rev_curr_hid;
        }
      } else {
        const VertexId vid = he_to(curr_hid);
        if (vid.valid()) {
          vertex_data(vid).halfedge = next_hid;
        }
      }

      if (rev_next_boundary) {
        remove_edge(he_edge(curr_hid));
      }
    }

    face_data(fid).halfedge = HalfedgeId{};
    n_faces_ = (n_faces_ == 0) ? 0 : (n_faces_ - 1);
  }

  FaceId new_face_by_halfedges(const std::span<const HalfedgeId> halfedges) {
    if (halfedges.empty()) {
      return FaceId{};
    }

    for (std::size_t i = 0; i < halfedges.size(); ++i) {
      const HalfedgeId ha_twin = halfedges[halfedges.size() - 1 - i];
      const HalfedgeId hb_twin = halfedges[halfedges.size() - 1 - ((i + 1) % halfedges.size())];
      const HalfedgeId ha = he_twin(ha_twin);
      const HalfedgeId hb = he_twin(hb_twin);

      const bool ha_boundary = he_is_boundary(ha);
      const bool hb_boundary = he_is_boundary(hb);

      if (ha_boundary && hb_boundary) {
        const VertexId vid = he_to(ha);
        const HalfedgeId vh = vid.valid() ? v_halfedge(vid) : HalfedgeId{};

        if (vid.valid()) {
          vertex_data(vid).halfedge = hb;
        }

        if (vh.valid()) {
          const HalfedgeId vh_prev = he_prev(vh);
          if (he_is_boundary(vh_prev) && he_is_boundary(vh)) {
            connect_halfedges(vh_prev, hb);
            connect_halfedges(ha, vh);
            continue;
          }
        }
        connect_halfedges(ha, hb);
      } else if (ha_boundary && !hb_boundary) {
        const HalfedgeId ha_next = he_next(hb_twin);
        connect_halfedges(ha, ha_next);
      } else if (!ha_boundary && hb_boundary) {
        const HalfedgeId hb_prev = he_prev(ha_twin);
        connect_halfedges(hb_prev, hb);
        const VertexId vid = he_to(ha);
        if (vid.valid()) {
          vertex_data(vid).halfedge = hb;
        }
      } else {
      }
    }

    const FaceId new_fid = new_face();
    for (std::size_t i = 0; i < halfedges.size(); ++i) {
      const HalfedgeId ha = halfedges[i];
      const HalfedgeId hb = halfedges[(i + 1) % halfedges.size()];
      connect_halfedges(ha, hb);
      halfedge_data(ha).face = new_fid;
    }
    face_data(new_fid).halfedge = halfedges[0];
    return new_fid;
  }

  FaceId new_face_by_halfedges(std::initializer_list<HalfedgeId> halfedges) {
    return new_face_by_halfedges(std::span<const HalfedgeId>{halfedges.begin(), halfedges.size()});
  }

  void split_face_into_triangles(const FaceId fid, std::span<const VertexId> triangles) noexcept {
    assert(triangles.size() % 3 == 0);
    std::unordered_map<std::pair<gpf::VertexId, gpf::VertexId>, HalfedgeId, detail::PairHash> edge_map;
    edge_map.reserve(std::ranges::distance(face(fid).halfedges()));

    for (const auto he : face(fid).halfedges()) {
      const VertexId va = he_from(he.id);
      const VertexId vb = he_to(he.id);
      edge_map.emplace((va.idx < vb.idx) ? std::pair{va, vb} : std::pair{vb, va}, he.id);
    }

    // Helper to get or create a halfedge from va to vb
    auto get_or_create_halfedge = [&](VertexId va, VertexId vb) -> HalfedgeId {
      auto key = (va.idx < vb.idx) ? std::pair{va, vb} : std::pair{vb, va};

      if (auto it = edge_map.find(key); it != edge_map.end()) {
        HalfedgeId hid = it->second;
        assert(he_to(hid) == vb);
        return hid;
      }

      // Create new edge
      HalfedgeId hid = new_edge();
      HalfedgeId twin_hid = he_twin(hid);

      halfedge_data(hid).vertex = vb;
      halfedge_data(twin_hid).vertex = va;

      auto& va_hid = vertex_data(va).halfedge;
      auto& vb_hid = vertex_data(vb).halfedge;
      if (!va_hid.valid()) {
        va_hid = hid;
      }
      if (!vb_hid.valid()) {
        vb_hid = twin_hid;
      }

      edge_map.emplace(std::move(key), twin_hid);
      return hid;
    };

    // Create triangular faces
    const std::size_t n_triangles = triangles.size() / 3;
    std::array<HalfedgeId, 3> hes;

    auto create_triangle = [this, &hes] (const auto new_fid) {
      auto ha = hes[0];
      auto hb = hes[1];
      auto hc = hes[2];
      halfedge_data(ha).face = new_fid;
      halfedge_data(hb).face = new_fid;
      halfedge_data(hc).face = new_fid;
      connect_halfedges(ha, hb);
      connect_halfedges(hb, hc);
      connect_halfedges(hc, ha);
      face_data(new_fid).halfedge = ha;
    };
    for (std::size_t i = 0; i < n_triangles; ++i) {
      const VertexId v0 = triangles[i * 3];
      const VertexId v1 = triangles[i * 3 + 1];
      const VertexId v2 = triangles[i * 3 + 2];

      hes[0] = get_or_create_halfedge(v0, v1);
      hes[1] = get_or_create_halfedge(v1, v2);
      hes[2] = get_or_create_halfedge(v2, v0);

      const auto new_fid = i == 0 ? fid : new_face();
      create_triangle(new_fid);
    }
  }

  auto he_prev_twin(HalfedgeId hid) const -> HalfedgeId {
    return halfedge(hid).prev().twin().id;
  }

  auto he_next_twin(HalfedgeId hid) const -> HalfedgeId {
    return halfedge(hid).next().twin().id;
  }

  auto he_twin_prev(HalfedgeId hid) const -> HalfedgeId {
    return halfedge(hid).twin().prev().id;
  }

  auto he_twin_next(HalfedgeId hid) const -> HalfedgeId {
      return halfedge(hid).twin().next().id;
  }

  auto he_to_to(HalfedgeId hid) const -> VertexId {
    return halfedge(hid).next().to().id;
  }

  void he_replace(HalfedgeId old_hid, HalfedgeId new_hid) {
    const VertexId va = he_from(old_hid);
    const HalfedgeId prev_hid = he_prev(old_hid);
    const HalfedgeId next_hid = he_next(old_hid);

    connect_halfedges(prev_hid, new_hid);
    connect_halfedges(new_hid, next_hid);

    const FaceId fid = he_face(old_hid);
    halfedge_data(new_hid).face = fid;

    if (va.valid() && v_halfedge(va) == old_hid) {
      vertex_data(va).halfedge = new_hid;
    }

    if (fid.valid()) {
      face_data(fid).halfedge = new_hid;
      halfedge_data(old_hid).face = FaceId{};
    }

    if (he_is_boundary(he_twin(old_hid))) {
      remove_edge(he_edge(old_hid));
    }
  }

  void flip(HalfedgeId hid) {
    const HalfedgeId bl_hid = he_next(hid);
    const HalfedgeId br_hid = he_next(bl_hid);

    const HalfedgeId twin_hid = he_twin(hid);
    const HalfedgeId tr_hid = he_next(twin_hid);
    const HalfedgeId tl_hid = he_next(tr_hid);

    const FaceId fid = he_face(hid);
    const FaceId twin_fid = he_face(twin_hid);

    const VertexId bottom_vid = he_to(bl_hid);
    const VertexId top_vid = he_to(tr_hid);
    const VertexId left_vid = he_to(hid);
    const VertexId right_vid = he_to(twin_hid);

    halfedge_data(tr_hid).face = fid;
    halfedge_data(bl_hid).face = twin_fid;

    halfedge_data(hid).vertex = bottom_vid;
    halfedge_data(twin_hid).vertex = top_vid;

    connect_halfedges(hid, br_hid);
    connect_halfedges(br_hid, tr_hid);
    connect_halfedges(tr_hid, hid);

    connect_halfedges(twin_hid, tl_hid);
    connect_halfedges(tl_hid, bl_hid);
    connect_halfedges(bl_hid, twin_hid);

    face_data(fid).halfedge = hid;
    face_data(twin_fid).halfedge = twin_hid;

    if (v_halfedge(left_vid) == twin_hid) {
      vertex_data(left_vid).halfedge = bl_hid;
    }
    if (v_halfedge(right_vid) == hid) {
      vertex_data(right_vid).halfedge = tr_hid;
    }
  }

  VertexId split_edge(EdgeId eid) {
    const HalfedgeId hid = e_halfedge(eid);
    const HalfedgeId twin_hid = he_twin(hid);
    const VertexId vb = he_to(hid);

    const VertexId new_v = new_vertices(1);
    const HalfedgeId new_hid = new_edge();
    const HalfedgeId new_twin_hid = he_twin(new_hid);

    if (vb.valid() && v_halfedge(vb) == twin_hid) {
      vertex_data(vb).halfedge = new_twin_hid;
    }
    vertex_data(new_v).halfedge = new_hid;

    const FaceId fid = he_face(hid);
    const FaceId twin_fid = he_face(twin_hid);

    halfedge_data(new_hid).face = fid;
    halfedge_data(new_twin_hid).face = twin_fid;

    halfedge_data(hid).vertex = new_v;
    halfedge_data(new_hid).vertex = vb;
    halfedge_data(new_twin_hid).vertex = new_v;

    const HalfedgeId prev_twin_hid = he_prev(twin_hid);
    const HalfedgeId next_hid = he_next(hid);

    connect_halfedges(prev_twin_hid, new_twin_hid);
    connect_halfedges(new_twin_hid, twin_hid);

    connect_halfedges(hid, new_hid);
    connect_halfedges(new_hid, next_hid);

    return new_v;
  }

  HalfedgeId split_face(FaceId fid, VertexId va, VertexId vb) {
    HalfedgeId left_last_hid{};
    HalfedgeId right_last_hid{};

    for (const auto he : face(fid).halfedges()) {
      const VertexId v = he.data().vertex;
      if (v == va) {
        left_last_hid = he.id;
      } else if (v == vb) {
        right_last_hid = he.id;
      }
    }

    assert(left_last_hid.valid());
    assert(right_last_hid.valid());

    const HalfedgeId left_first_hid = he_next(right_last_hid);
    const HalfedgeId right_first_hid = he_next(left_last_hid);

    const HalfedgeId first_he = new_edge();
    const HalfedgeId second_he = he_twin(first_he);

    halfedge_data(first_he).vertex = va;
    halfedge_data(second_he).vertex = vb;

    connect_halfedges(right_last_hid, first_he);
    connect_halfedges(first_he, right_first_hid);
    connect_halfedges(left_last_hid, second_he);
    connect_halfedges(second_he, left_first_hid);

    const FaceId new_f = new_face();
    halfedge_data(first_he).face = fid;
    {
      HalfedgeId hid = second_he;
      do {
        halfedge_data(hid).face = new_f;
        hid = he_next(hid);
      } while (hid != second_he);
    }

    face_data(fid).halfedge = first_he;
    face_data(new_f).halfedge = second_he;
    return second_he;
  }

 private:
  std::vector<VertexData> vertices_;
  std::vector<HalfedgeData> halfedges_;
  std::vector<FaceData> faces_;
  std::vector<EdgeProp> edges_;

  std::size_t n_vertices_{0ull};
  std::size_t n_halfedges_{0ull};
  std::size_t n_faces_{0ull};

  void v_min_reserve(VertexId vid) {
    if (!vid.valid()) {
      return;
    }
    const std::size_t len = vid.idx + 1;
    if (vertices_.size() < len) {
      vertices_.resize(len);
    }
  }

  [[nodiscard]] bool he_is_valid(HalfedgeId hid) const { return halfedge_data(hid).vertex.valid(); }

  [[nodiscard]] bool e_is_valid(EdgeId eid) const {
    const std::size_t idx = eid.idx << 1;
    return he_is_valid(HalfedgeId{idx}) || he_is_valid(HalfedgeId{idx + 1});
  }

  void recount_n_vertices() {
    n_vertices_ = std::ranges::count_if(vertices_, [](const VertexData& v) { return v.valid(); });
  }

  template <class PolygonRange>
  void build(PolygonRange&& polygons) {
    std::unordered_map<std::pair<std::size_t, std::size_t>, HalfedgeId, detail::PairHash> edge_map;

    for (const auto& polygon : polygons) {
      std::vector<std::size_t> verts;
      for (const auto v_raw : polygon) {
        verts.push_back(static_cast<std::size_t>(v_raw));
      }

      if (verts.empty()) {
        continue;
      }

      const FaceId fid = new_face();
      HalfedgeId first_hid{};
      HalfedgeId prev_hid{};

      auto add_edge = [&](const std::size_t a, const std::size_t b) {
        const auto key = (a < b) ? std::pair{a, b} : std::pair{b, a};
        const VertexId va{a};
        const VertexId vb{b};

        if (va.valid()) {
          v_min_reserve(va);
        }
        if (vb.valid()) {
          v_min_reserve(vb);
        }

        HalfedgeId hid{};
        if (auto it = edge_map.find(key); it != edge_map.end()) {
          hid = he_twin(it->second);
          it->second = hid;
        } else {
          hid = new_edge();
          edge_map.emplace(key, hid);
        }

        if (va.valid()) {
          vertex_data(va).halfedge = hid;
        }

        halfedge_data(hid).vertex = vb;
        halfedge_data(he_twin(hid)).vertex = va;
        halfedge_data(hid).face = fid;

        if (!first_hid.valid()) {
          face_data(fid).halfedge = hid;
          first_hid = hid;
        } else {
          connect_halfedges(prev_hid, hid);
        }
        prev_hid = hid;
      };

      for (std::size_t i = 0; i < verts.size(); ++i) {
        const std::size_t a = verts[i];
        const std::size_t b = verts[(i + 1) % verts.size()];
        add_edge(a, b);
      }

      connect_halfedges(prev_hid, first_hid);
    }

    recount_n_vertices();

    std::vector<bool> edge_visited(edges_.size(), false);
    for (const auto& [_, twin_hid] : edge_map) {
      if ((twin_hid.idx & 1) != 0) {
        continue;
      }

      const EdgeId eid = he_edge(twin_hid);
      if (edge_visited[eid.idx]) {
        continue;
      }
      edge_visited[eid.idx] = true;

      const HalfedgeId first_hid = he_twin(twin_hid);
      assert(he_is_boundary(first_hid));

      HalfedgeId curr_hid = first_hid;
      while (true) {
        HalfedgeId prev_hid = curr_hid;
        while (true) {
          prev_hid = he_next(he_twin(prev_hid));
          const HalfedgeId twin_prev_hid = he_twin(prev_hid);
          if (he_is_boundary(twin_prev_hid)) {
            prev_hid = twin_prev_hid;
            break;
          }
        }

        const VertexId vid = he_to(prev_hid);
        if (vid.valid()) {
          vertex_data(vid).halfedge = curr_hid;
        }

        edge_visited[he_edge(prev_hid).idx] = true;
        connect_halfedges(prev_hid, curr_hid);
        curr_hid = prev_hid;
        if (curr_hid == first_hid) {
          break;
        }
      }
    }
  }
};

}  // namespace gpf
