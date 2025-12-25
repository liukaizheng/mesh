#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <ranges>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gpf/detail.hpp"
#include "gpf/handles.hpp"
#include "gpf/ids.hpp"

namespace gpf {

template <class VertexProp, class HalfedgeProp, class EdgeProp, class FaceProp>
class SurfaceMesh {
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

    EdgeId edge{};
    HalfedgeId sibling{};
    HalfedgeId incoming_next{};

    HalfedgeProp property{};
    [[nodiscard]] bool valid() const { return vertex.valid(); }
  };

  struct FaceData {
    HalfedgeId halfedge{};
    FaceProp property{};
    [[nodiscard]] bool valid() const { return halfedge.valid(); }
  };

  struct EdgeData {
    HalfedgeId halfedge{};
    EdgeProp property{};
    [[nodiscard]] bool valid() const { return halfedge.valid(); }
  };

  using Self = SurfaceMesh<VertexProp, HalfedgeProp, EdgeProp, FaceProp>;
  using Vertex = VertexHandle<Self, false>;
  using Halfedge = HalfedgeHandle<Self, false>;
  using Face = FaceHandle<Self, false>;
  using Edge = EdgeHandle<Self, false>;
  using ConstVertex = VertexHandle<Self, true>;
  using ConstHalfedge = HalfedgeHandle<Self, true>;
  using ConstFace = FaceHandle<Self, true>;
  using ConstEdge = EdgeHandle<Self, true>;

  SurfaceMesh() = default;

  template <class PolygonRange>
  static Self new_in(PolygonRange&& polygons) {
    Self mesh;
    mesh.build(std::forward<PolygonRange>(polygons));
    return mesh;
  }

  [[nodiscard]] std::size_t n_vertices() const { return n_vertices_; }
  [[nodiscard]] std::size_t n_halfedges() const { return n_halfedges_; }
  [[nodiscard]] std::size_t n_edges() const { return n_edges_; }
  [[nodiscard]] std::size_t n_faces() const { return n_faces_; }

  [[nodiscard]] std::size_t n_vertices_capacity() const { return vertices_.size(); }
  [[nodiscard]] std::size_t n_halfedges_capacity() const { return halfedges_.size(); }
  [[nodiscard]] std::size_t n_edges_capacity() const { return edges_.size(); }
  [[nodiscard]] std::size_t n_faces_capacity() const { return faces_.size(); }

  [[nodiscard]] const VertexData& vertex_data(VertexId vid) const { return vertices_[vid.idx]; }
  [[nodiscard]] VertexData& vertex_data(VertexId vid) { return vertices_[vid.idx]; }

  [[nodiscard]] const HalfedgeData& halfedge_data(HalfedgeId hid) const { return halfedges_[hid.idx]; }
  [[nodiscard]] HalfedgeData& halfedge_data(HalfedgeId hid) { return halfedges_[hid.idx]; }

  [[nodiscard]] const FaceData& face_data(FaceId fid) const { return faces_[fid.idx]; }
  [[nodiscard]] FaceData& face_data(FaceId fid) { return faces_[fid.idx]; }

  [[nodiscard]] const EdgeData& edge_data(EdgeId eid) const { return edges_[eid.idx]; }
  [[nodiscard]] EdgeData& edge_data(EdgeId eid) { return edges_[eid.idx]; }

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
           std::views::filter([this](const std::size_t i) { return edges_[i].valid(); }) |
           std::views::transform([this](const std::size_t i) { return Edge{EdgeId{i}, this}; });
  }

  [[nodiscard]] auto edges() const {
    return std::views::iota(std::size_t{0}, edges_.size()) |
           std::views::filter([this](const std::size_t i) { return edges_[i].valid(); }) |
           std::views::transform([this](const std::size_t i) { return ConstEdge{EdgeId{i}, this}; });
  }

  [[nodiscard]] VertexId he_from(HalfedgeId hid) const { return halfedges_[halfedges_[hid.idx].prev.idx].vertex; }
  [[nodiscard]] VertexId he_to(HalfedgeId hid) const { return halfedges_[hid.idx].vertex; }
  [[nodiscard]] std::array<VertexId, 2> he_vertices(HalfedgeId hid) const {
    return {he_from(hid), he_to(hid)};
  }
  [[nodiscard]] HalfedgeId he_prev(HalfedgeId hid) const { return halfedges_[hid.idx].prev; }
  [[nodiscard]] HalfedgeId he_next(HalfedgeId hid) const { return halfedges_[hid.idx].next; }
  [[nodiscard]] FaceId he_face(HalfedgeId hid) const { return halfedges_[hid.idx].face; }

  [[nodiscard]] EdgeId he_edge(HalfedgeId hid) const { return halfedges_[hid.idx].edge; }
  [[nodiscard]] HalfedgeId he_sibling(HalfedgeId hid) const { return halfedges_[hid.idx].sibling; }
  [[nodiscard]] HalfedgeId he_incoming_next(HalfedgeId hid) const { return halfedges_[hid.idx].incoming_next; }

  [[nodiscard]] HalfedgeId e_halfedge(EdgeId eid) const { return edges_[eid.idx].halfedge; }
  [[nodiscard]] std::array<VertexId, 2> e_vertices(EdgeId eid) const {
    const HalfedgeId hid = e_halfedge(eid);
    return he_vertices(hid);
  }

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

  HalfedgeId new_halfedges(std::size_t n) {
    const HalfedgeId ret{halfedges_.size()};
    halfedges_.resize(ret.idx + n);
    n_halfedges_ += n;
    return ret;
  }

  FaceId new_faces(std::size_t n) {
    const FaceId ret{faces_.size()};
    faces_.resize(ret.idx + n);
    n_faces_ += n;
    return ret;
  }

  EdgeId new_edges(std::size_t n) {
    const EdgeId ret{edges_.size()};
    edges_.resize(ret.idx + n);
    n_edges_ += n;
    return ret;
  }

  VertexId split_edge(EdgeId eid) {
    const HalfedgeId edge_hid = e_halfedge(eid);
    const VertexId vb = he_to(edge_hid);

    std::vector<HalfedgeId> old_halfedges;
    {
      HalfedgeId curr = edge_hid;
      do {
        old_halfedges.push_back(curr);
        curr = he_sibling(curr);
      } while (curr != edge_hid);
    }

    const std::size_t n_edge_halfedges = old_halfedges.size();
    const HalfedgeId start_hid = new_halfedges(n_edge_halfedges);
    const EdgeId new_eid = new_edges(1);
    const VertexId new_vid = new_vertices(1);
    vertex_data(new_vid).halfedge = edge_hid;

    HalfedgeId old_first{};
    HalfedgeId new_first{};
    HalfedgeId old_prev{};
    HalfedgeId new_prev{};

    for (std::size_t i = 0; i < n_edge_halfedges; ++i) {
      const HalfedgeId old_hid = old_halfedges[i];
      const HalfedgeId new_hid = start_hid + i;

      auto& old_he = halfedge_data(old_hid);
      auto& new_he = halfedge_data(new_hid);

      new_he.vertex = new_vid;
      new_he.face = old_he.face;

      const VertexId from_vid = he_from(old_hid);
      if (from_vid.valid() && vertex_data(from_vid).halfedge == old_hid) {
        vertex_data(from_vid).halfedge = new_hid;
      }

      const HalfedgeId prev_hid = old_he.prev;
      connect_halfedges(prev_hid, new_hid);
      connect_halfedges(new_hid, old_hid);

      HalfedgeId old_edge_he = old_hid;
      HalfedgeId new_edge_he = new_hid;
      if (from_vid == vb) {
        old_he.edge = new_eid;
        new_he.edge = eid;
        old_edge_he = new_hid;
        new_edge_he = old_hid;
      } else {
        new_he.edge = new_eid;
      }

      if (!old_first.valid()) {
        old_first = old_edge_he;
        new_first = new_edge_he;
      } else {
        halfedge_data(old_prev).sibling = old_edge_he;
        halfedge_data(new_prev).sibling = new_edge_he;

        if (halfedge_data(old_prev).vertex == new_vid) {
          halfedge_data(old_prev).incoming_next = new_hid;
        } else {
          halfedge_data(new_prev).incoming_next = new_hid;
        }
      }

      old_prev = old_edge_he;
      new_prev = new_edge_he;
    }

    edge_data(eid).halfedge = old_first;
    edge_data(new_eid).halfedge = new_first;

    if (halfedge_data(old_prev).vertex == new_vid) {
      halfedge_data(old_prev).incoming_next = start_hid;
    } else {
      halfedge_data(new_prev).incoming_next = start_hid;
    }

    halfedge_data(old_prev).sibling = old_first;
    halfedge_data(new_prev).sibling = new_first;

    return new_vid;
  }

  HalfedgeId split_face(FaceId fid, VertexId va, VertexId vb) {
    const HalfedgeId new_start_hid = new_halfedges(2);
    const EdgeId new_eid = new_edges(1);
    const FaceId new_fid = new_faces(1);

    HalfedgeId left_last_hid{};
    HalfedgeId right_last_hid{};
    for (const auto he : face(fid).halfedges()) {
      if (he.data().vertex == va) {
        left_last_hid = he.id;
      } else if (he.data().vertex == vb) {
        right_last_hid = he.id;
      }
    }

    assert(left_last_hid.valid());
    assert(right_last_hid.valid());

    const HalfedgeId left_first_hid = he_next(right_last_hid);
    const HalfedgeId right_first_hid = he_next(left_last_hid);

    const HalfedgeId first_hid = new_start_hid;
    const HalfedgeId second_hid = new_start_hid + 1;

    halfedge_data(first_hid).sibling = second_hid;
    halfedge_data(second_hid).sibling = first_hid;

    connect_halfedges(right_last_hid, first_hid);
    connect_halfedges(first_hid, right_first_hid);
    connect_halfedges(left_last_hid, second_hid);
    connect_halfedges(second_hid, left_first_hid);

    auto& first_he = halfedge_data(first_hid);
    first_he.vertex = va;
    first_he.edge = new_eid;
    first_he.face = fid;

    auto& second_he = halfedge_data(second_hid);
    second_he.vertex = vb;
    second_he.edge = new_eid;
    second_he.face = new_fid;

    edge_data(new_eid).halfedge = second_hid;

    if (he_to(left_last_hid).valid()) {
      const HalfedgeId next = he_incoming_next(left_last_hid);
      halfedge_data(left_last_hid).incoming_next = first_hid;
      halfedge_data(first_hid).incoming_next = next;
    }
    if (he_to(right_last_hid).valid()) {
      const HalfedgeId next = he_incoming_next(right_last_hid);
      halfedge_data(right_last_hid).incoming_next = second_hid;
      halfedge_data(second_hid).incoming_next = next;
    }

    {
      HalfedgeId curr = left_first_hid;
      while (curr != second_hid) {
        halfedge_data(curr).face = new_fid;
        curr = he_next(curr);
      }
    }

    face_data(fid).halfedge = first_hid;
    face_data(new_fid).halfedge = second_hid;
    return second_hid;
  }

  FaceId add_face_by_halfedges(const std::vector<HalfedgeId>& halfedges, bool build_siblings) {
    if (halfedges.empty()) {
      return FaceId{};
    }

    const HalfedgeId first_new_hid = new_halfedges(halfedges.size());
    const FaceId new_fid = new_faces(1);

    HalfedgeId prev_new_hid{};
    for (std::size_t i = 0; i < halfedges.size(); ++i) {
      const HalfedgeId old_hid = halfedges[i];
      const HalfedgeId new_hid = first_new_hid + i;

      const auto& old_he = halfedge_data(old_hid);
      auto& new_he = halfedge_data(new_hid);
      new_he.vertex = old_he.vertex;
      new_he.edge = old_he.edge;
      new_he.face = new_fid;

      if (build_siblings) {
        const HalfedgeId old_next = he_sibling(old_hid);
        halfedge_data(old_hid).sibling = new_hid;
        new_he.sibling = old_next;
      }
      if (old_he.vertex.valid()) {
        const HalfedgeId old_next = he_incoming_next(old_hid);
        halfedge_data(old_hid).incoming_next = new_hid;
        new_he.incoming_next = old_next;
      }

      if (i != 0) {
        connect_halfedges(prev_new_hid, new_hid);
      }
      prev_new_hid = new_hid;
    }
    connect_halfedges(prev_new_hid, first_new_hid);
    face_data(new_fid).halfedge = first_new_hid;
    return new_fid;
  }

  void set_he_sibling(HalfedgeId hid, HalfedgeId sid) { halfedge_data(hid).sibling = sid; }

 private:
  std::vector<VertexData> vertices_;
  std::vector<HalfedgeData> halfedges_;
  std::vector<FaceData> faces_;
  std::vector<EdgeData> edges_;

  std::size_t n_vertices_{0};
  std::size_t n_halfedges_{0};
  std::size_t n_faces_{0};
  std::size_t n_edges_{0};

  void v_min_reserve(VertexId vid) {
    if (!vid.valid()) {
      return;
    }
    const std::size_t len = vid.idx + 1;
    if (vertices_.size() < len) {
      vertices_.resize(len);
    }
  }

  void recount_n_vertices() {
    n_vertices_ = std::ranges::count_if(vertices_, [](const VertexData& v) { return v.valid(); });
  }

  std::pair<std::vector<HalfedgeId>, std::vector<std::size_t>> vertex_cycle() const {
    std::vector<std::size_t> v_degree(vertices_.size(), 0);
    for (const auto& he : halfedges_) {
      if (he.vertex.valid()) {
        v_degree[he.vertex.idx] += 1;
      }
    }

    std::vector<std::size_t> vertex_separators;
    vertex_separators.reserve(vertices_.size() + 1);
    vertex_separators.push_back(0);
    std::size_t sum = 0;
    for (const auto degree : v_degree) {
      sum += degree;
      vertex_separators.push_back(sum);
    }

    std::vector<std::size_t> he_positions = vertex_separators;
    std::vector<HalfedgeId> vertex_halfedges(halfedges_.size(), HalfedgeId{0});

    for (std::size_t hid = 0; hid < halfedges_.size(); ++hid) {
      const auto& he = halfedges_[hid];
      if (!he.vertex.valid()) {
        continue;
      }
      const std::size_t pos = he_positions[he.vertex.idx];
      vertex_halfedges[pos] = HalfedgeId{hid};
      he_positions[he.vertex.idx] += 1;
    }

    return {std::move(vertex_halfedges), std::move(vertex_separators)};
  }

  template <class PolygonRange>
  void build(PolygonRange&& polygons) {
    for (const auto& polygon : polygons) {
      bool first = true;

      FaceId fid{};
      HalfedgeId first_hid{};
      HalfedgeId prev_hid{};
      VertexId prev_vid{};

      for (const auto v_raw : polygon) {
        const VertexId vid{static_cast<std::size_t>(v_raw)};
        if (vid.valid()) {
          v_min_reserve(vid);
        }

        const HalfedgeId hid = new_halfedges(1);
        auto& he = halfedge_data(hid);
        he.vertex = vid;
        if (first) {
          fid = new_faces(1);
          he.face = fid;
          face_data(fid).halfedge = hid;
          first_hid = hid;
          first = false;
        } else {
          he.face = fid;
          if (prev_vid.valid()) {
            vertex_data(prev_vid).halfedge = hid;
          }
          connect_halfedges(prev_hid, hid);
        }
        prev_vid = vid;
        prev_hid = hid;
      }

      if (first) {
        continue;
      }
      if (prev_vid.valid()) {
        vertex_data(prev_vid).halfedge = first_hid;
      }
      connect_halfedges(prev_hid, first_hid);
    }

    recount_n_vertices();

    std::unordered_map<std::pair<VertexId, VertexId>, HalfedgeId, detail::PairHash> edge_history;

    for (std::size_t hid = 0; hid < halfedges_.size(); ++hid) {
      const HalfedgeId he_id{hid};
      const auto [va, vb] = he_vertices(he_id);
      const auto key = detail::ordered_pair(va, vb);

      if (auto it = edge_history.find(key); it != edge_history.end()) {
        const HalfedgeId prev_seen = it->second;
        const EdgeId eid = he_edge(prev_seen);
        auto& he = halfedge_data(he_id);
        he.sibling = prev_seen;
        he.edge = eid;
        it->second = he_id;
      } else {
        const EdgeId new_eid = new_edges(1);
        auto& he = halfedge_data(he_id);
        he.edge = new_eid;
        he.sibling = HalfedgeId{};
        edge_data(new_eid).halfedge = he_id;
        edge_history.emplace(key, he_id);
      }
    }

    for (const auto& [_, last_hid] : edge_history) {
      HalfedgeId curr = last_hid;
      while (halfedge_data(curr).sibling.valid()) {
        curr = halfedge_data(curr).sibling;
      }
      halfedge_data(curr).sibling = last_hid;
    }

    const auto [v_in_halfedges, v_in_separators] = vertex_cycle();
    for (std::size_t vid = 0; vid + 1 < v_in_separators.size(); ++vid) {
      const std::size_t start = v_in_separators[vid];
      const std::size_t end = v_in_separators[vid + 1];
      if (end <= start) {
        continue;
      }

      for (std::size_t i = start; i < end; ++i) {
        const HalfedgeId ha = v_in_halfedges[i];
        const HalfedgeId hb = v_in_halfedges[(i + 1 == end) ? start : (i + 1)];
        halfedge_data(ha).incoming_next = hb;
      }
    }
  }
};

}  // namespace gpf
