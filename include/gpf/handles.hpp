#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <type_traits>

#include "gpf/ids.hpp"

namespace gpf {

template <class Mesh, bool Const>
using mesh_ptr_t = std::conditional_t<Const, const Mesh*, Mesh*>;

template <class Mesh, bool Const>
class VertexHandle;
template <class Mesh, bool Const>
class HalfedgeHandle;
template <class Mesh, bool Const>
class FaceHandle;
template <class Mesh, bool Const>
class EdgeHandle;

template <class Mesh, bool Const, HalfedgeId (Mesh::*Next)(HalfedgeId) const>
class HalfedgeCycleRange;

template <class Mesh, bool Const>
class VertexNeighborEdgesRange;
template <class Mesh, bool Const>
class VertexNeighborsRange;

template <class Mesh, bool Const>
class VertexHandle {
 public:
  using mesh_ptr = mesh_ptr_t<Mesh, Const>;
  using data_ref = std::conditional_t<Const, const typename Mesh::VertexData&, typename Mesh::VertexData&>;

  VertexId id{};
  mesh_ptr mesh{};

  VertexHandle() = default;
  VertexHandle(VertexId id_, mesh_ptr mesh_) : id(id_), mesh(mesh_) {}

  [[nodiscard]] data_ref data() const { return mesh->vertex_data(id); }

  [[nodiscard]] HalfedgeHandle<Mesh, Const> halfedge() const;
  [[nodiscard]] auto incoming_halfedges() const
      -> HalfedgeCycleRange<Mesh, Const, &Mesh::he_incoming_next>;
  [[nodiscard]] auto outgoing_halfedges() const;

  [[nodiscard]] auto edges() const -> VertexNeighborEdgesRange<Mesh, Const>;
  [[nodiscard]] auto vertices() const -> VertexNeighborsRange<Mesh, Const>;
};

template <class Mesh, bool Const>
class FaceHandle {
 public:
  using mesh_ptr = mesh_ptr_t<Mesh, Const>;
  using data_ref = std::conditional_t<Const, const typename Mesh::FaceData&, typename Mesh::FaceData&>;

  FaceId id{};
  mesh_ptr mesh{};

  FaceHandle() = default;
  FaceHandle(FaceId id_, mesh_ptr mesh_) : id(id_), mesh(mesh_) {}

  [[nodiscard]] data_ref data() const { return mesh->face_data(id); }

  [[nodiscard]] HalfedgeHandle<Mesh, Const> halfedge() const;

  [[nodiscard]] auto halfedges() const -> HalfedgeCycleRange<Mesh, Const, &Mesh::he_next>;
  [[nodiscard]] auto halfedges_reverse() const -> HalfedgeCycleRange<Mesh, Const, &Mesh::he_prev>;
};

template <class Mesh, bool Const>
class HalfedgeHandle {
 public:
  using mesh_ptr = mesh_ptr_t<Mesh, Const>;
  using data_ref =
      std::conditional_t<Const, const typename Mesh::HalfedgeData&, typename Mesh::HalfedgeData&>;

  HalfedgeId id{};
  mesh_ptr mesh{};

  HalfedgeHandle() = default;
  HalfedgeHandle(HalfedgeId id_, mesh_ptr mesh_) : id(id_), mesh(mesh_) {}

  [[nodiscard]] data_ref data() const { return mesh->halfedge_data(id); }

  [[nodiscard]] VertexHandle<Mesh, Const> from() const;
  [[nodiscard]] VertexHandle<Mesh, Const> to() const;
  [[nodiscard]] HalfedgeHandle<Mesh, Const> next() const;
  [[nodiscard]] HalfedgeHandle<Mesh, Const> prev() const;
  [[nodiscard]] FaceHandle<Mesh, Const> face() const;

  [[nodiscard]] EdgeHandle<Mesh, Const> edge() const;
  [[nodiscard]] HalfedgeHandle<Mesh, Const> sibling() const;
  [[nodiscard]] HalfedgeHandle<Mesh, Const> incoming_next() const;
  [[nodiscard]] HalfedgeHandle<Mesh, Const> twin() const;
};

template <class Mesh, bool Const>
class EdgeHandle {
 public:
  using mesh_ptr = mesh_ptr_t<Mesh, Const>;
  using data_ref = std::conditional_t<Const, const typename Mesh::EdgeData&, typename Mesh::EdgeData&>;

  EdgeId id{};
  mesh_ptr mesh{};

  EdgeHandle() = default;
  EdgeHandle(EdgeId id_, mesh_ptr mesh_) : id(id_), mesh(mesh_) {}

  [[nodiscard]] data_ref data() const { return mesh->edge_data(id); }

  [[nodiscard]] HalfedgeHandle<Mesh, Const> halfedge() const;
  [[nodiscard]] auto halfedges() const -> HalfedgeCycleRange<Mesh, Const, &Mesh::he_sibling>;
  [[nodiscard]] std::array<VertexHandle<Mesh, Const>, 2> vertices() const;
};

template <class Mesh, bool Const, HalfedgeId (Mesh::*Next)(HalfedgeId) const>
class HalfedgeCycleRange {
 public:
  using mesh_ptr = mesh_ptr_t<Mesh, Const>;
  using value_type = HalfedgeHandle<Mesh, Const>;

  HalfedgeCycleRange(mesh_ptr mesh, HalfedgeId start) : mesh_(mesh), start_(start) {}

  class iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = HalfedgeHandle<Mesh, Const>;

    iterator() = default;
    iterator(mesh_ptr mesh, HalfedgeId first, bool done)
        : mesh_(mesh), first_(first), current_(first), first_step_(true), done_(done) {}

    [[nodiscard]] value_type operator*() const { return value_type{current_, mesh_}; }

    iterator& operator++() {
      if (done_) {
        return *this;
      }

      first_step_ = false;
      current_ = (mesh_->*Next)(current_);
      if (current_ == first_) {
        done_ = true;
      }
      return *this;
    }

    void operator++(int) { (void)operator++(); }

    [[nodiscard]] bool operator==(std::default_sentinel_t) const { return done_; }

   private:
    mesh_ptr mesh_{};
    HalfedgeId first_{};
    HalfedgeId current_{};
    bool first_step_{true};
    bool done_{true};
  };

  [[nodiscard]] iterator begin() const {
    if (!start_.valid()) {
      return iterator{mesh_, start_, true};
    }
    return iterator{mesh_, start_, false};
  }

  [[nodiscard]] std::default_sentinel_t end() const { return {}; }

 private:
  mesh_ptr mesh_{};
  HalfedgeId start_{};
};

template <class Mesh, bool Const>
class VertexNeighborEdgesRange {
 public:
  using mesh_ptr = mesh_ptr_t<Mesh, Const>;
  using value_type = EdgeHandle<Mesh, Const>;

  VertexNeighborEdgesRange(mesh_ptr mesh, VertexId vid) : mesh_(mesh), vid_(vid) {}

  class iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = EdgeHandle<Mesh, Const>;

    iterator() = default;
    iterator(mesh_ptr mesh, VertexId vid, bool end) : mesh_(mesh), done_(end) {
      if (done_) {
        return;
      }
      const HalfedgeId out_hid = mesh_->vertex_data(vid).halfedge;
      if (!out_hid.valid()) {
        done_ = true;
        return;
      }

      first_hid_ = mesh_->he_prev(out_hid);
      hid_ = first_hid_;
      visited_ = false;
      first_ = true;
      advance();
    }

    [[nodiscard]] value_type operator*() const {
      assert(curr_eid_.valid());
      return value_type{curr_eid_, mesh_};
    }

    iterator& operator++() {
      advance();
      return *this;
    }

    void operator++(int) { (void)operator++(); }

    [[nodiscard]] bool operator==(std::default_sentinel_t) const { return done_; }

   private:
    [[nodiscard]] bool valid() const { return first_ || visited_ || hid_ != first_hid_; }

    void next_halfedge() {
      first_ = false;
      if (visited_) {
        visited_ = false;
        hid_ = mesh_->he_incoming_next(hid_);
      } else {
        visited_ = true;
      }
    }

    void advance() {
      while (true) {
        if (!valid()) {
          done_ = true;
          curr_eid_ = EdgeId{};
          return;
        }

        if (visited_) {
          const HalfedgeId out_hid = mesh_->he_next(hid_);
          const EdgeId eid = mesh_->he_edge(out_hid);
          if (mesh_->e_halfedge(eid) == out_hid) {
            curr_eid_ = eid;
            next_halfedge();
            return;
          }
        } else {
          const EdgeId eid = mesh_->he_edge(hid_);
          if (mesh_->e_halfedge(eid) == hid_) {
            curr_eid_ = eid;
            next_halfedge();
            return;
          }
        }

        next_halfedge();
      }
    }

    mesh_ptr mesh_{};
    HalfedgeId first_hid_{};
    HalfedgeId hid_{};
    bool visited_{false};
    bool first_{true};
    bool done_{true};
    EdgeId curr_eid_{};
  };

  [[nodiscard]] iterator begin() const { return iterator{mesh_, vid_, false}; }
  [[nodiscard]] std::default_sentinel_t end() const { return {}; }

 private:
  mesh_ptr mesh_{};
  VertexId vid_{};
};

template <class Mesh, bool Const>
class VertexNeighborsRange {
 public:
  using mesh_ptr = mesh_ptr_t<Mesh, Const>;
  using value_type = VertexHandle<Mesh, Const>;

  VertexNeighborsRange(mesh_ptr mesh, VertexId vid) : mesh_(mesh), vid_(vid) {}

  class iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = VertexHandle<Mesh, Const>;

    iterator() = default;
    iterator(mesh_ptr mesh, VertexId vid, bool end) : mesh_(mesh), done_(end) {
      if (done_) {
        return;
      }
      const HalfedgeId out_hid = mesh_->vertex_data(vid).halfedge;
      if (!out_hid.valid()) {
        done_ = true;
        return;
      }

      first_hid_ = mesh_->he_prev(out_hid);
      hid_ = first_hid_;
      visited_ = false;
      first_ = true;
      advance();
    }

    [[nodiscard]] value_type operator*() const {
      assert(curr_vid_.valid());
      return value_type{curr_vid_, mesh_};
    }

    iterator& operator++() {
      advance();
      return *this;
    }

    void operator++(int) { (void)operator++(); }

    [[nodiscard]] bool operator==(std::default_sentinel_t) const { return done_; }

   private:
    [[nodiscard]] bool valid() const { return first_ || visited_ || hid_ != first_hid_; }

    void next_halfedge() {
      first_ = false;
      if (visited_) {
        visited_ = false;
        hid_ = mesh_->he_incoming_next(hid_);
      } else {
        visited_ = true;
      }
    }

    void advance() {
      while (true) {
        if (!valid()) {
          done_ = true;
          curr_vid_ = VertexId{};
          return;
        }

        if (visited_) {
          const HalfedgeId out_hid = mesh_->he_next(hid_);
          const EdgeId eid = mesh_->he_edge(out_hid);
          if (mesh_->e_halfedge(eid) == out_hid) {
            curr_vid_ = mesh_->he_to(out_hid);
            next_halfedge();
            return;
          }
        } else {
          const EdgeId eid = mesh_->he_edge(hid_);
          if (mesh_->e_halfedge(eid) == hid_) {
            curr_vid_ = mesh_->he_from(hid_);
            next_halfedge();
            return;
          }
        }

        next_halfedge();
      }
    }

    mesh_ptr mesh_{};
    HalfedgeId first_hid_{};
    HalfedgeId hid_{};
    bool visited_{false};
    bool first_{true};
    bool done_{true};
    VertexId curr_vid_{};
  };

  [[nodiscard]] iterator begin() const { return iterator{mesh_, vid_, false}; }
  [[nodiscard]] std::default_sentinel_t end() const { return {}; }

 private:
  mesh_ptr mesh_{};
  VertexId vid_{};
};

template <class Mesh, bool Const>
HalfedgeHandle<Mesh, Const> VertexHandle<Mesh, Const>::halfedge() const {
  return HalfedgeHandle<Mesh, Const>{data().halfedge, mesh};
}

template <class Mesh, bool Const>
auto VertexHandle<Mesh, Const>::incoming_halfedges() const
    -> HalfedgeCycleRange<Mesh, Const, &Mesh::he_incoming_next> {
  const HalfedgeId out_hid = data().halfedge;
  const HalfedgeId start = out_hid.valid() ? mesh->he_prev(out_hid) : HalfedgeId{};
  return HalfedgeCycleRange<Mesh, Const, &Mesh::he_incoming_next>{mesh, start};
}

template <class Mesh, bool Const>
auto VertexHandle<Mesh, Const>::outgoing_halfedges() const {
  return incoming_halfedges() |
         std::views::transform([](const HalfedgeHandle<Mesh, Const> he) { return he.next(); });
}

template <class Mesh, bool Const>
auto VertexHandle<Mesh, Const>::edges() const -> VertexNeighborEdgesRange<Mesh, Const> {
  return VertexNeighborEdgesRange<Mesh, Const>{mesh, id};
}

template <class Mesh, bool Const>
auto VertexHandle<Mesh, Const>::vertices() const -> VertexNeighborsRange<Mesh, Const> {
  return VertexNeighborsRange<Mesh, Const>{mesh, id};
}

template <class Mesh, bool Const>
HalfedgeHandle<Mesh, Const> FaceHandle<Mesh, Const>::halfedge() const {
  return HalfedgeHandle<Mesh, Const>{data().halfedge, mesh};
}

template <class Mesh, bool Const>
auto FaceHandle<Mesh, Const>::halfedges() const -> HalfedgeCycleRange<Mesh, Const, &Mesh::he_next> {
  return HalfedgeCycleRange<Mesh, Const, &Mesh::he_next>{mesh, data().halfedge};
}

template <class Mesh, bool Const>
auto FaceHandle<Mesh, Const>::halfedges_reverse() const
    -> HalfedgeCycleRange<Mesh, Const, &Mesh::he_prev> {
  return HalfedgeCycleRange<Mesh, Const, &Mesh::he_prev>{mesh, data().halfedge};
}

template <class Mesh, bool Const>
VertexHandle<Mesh, Const> HalfedgeHandle<Mesh, Const>::from() const {
  return VertexHandle<Mesh, Const>{mesh->he_from(id), mesh};
}

template <class Mesh, bool Const>
VertexHandle<Mesh, Const> HalfedgeHandle<Mesh, Const>::to() const {
  return VertexHandle<Mesh, Const>{data().vertex, mesh};
}

template <class Mesh, bool Const>
HalfedgeHandle<Mesh, Const> HalfedgeHandle<Mesh, Const>::next() const {
  return HalfedgeHandle<Mesh, Const>{data().next, mesh};
}

template <class Mesh, bool Const>
HalfedgeHandle<Mesh, Const> HalfedgeHandle<Mesh, Const>::prev() const {
  return HalfedgeHandle<Mesh, Const>{data().prev, mesh};
}

template <class Mesh, bool Const>
FaceHandle<Mesh, Const> HalfedgeHandle<Mesh, Const>::face() const {
  return FaceHandle<Mesh, Const>{data().face, mesh};
}

template <class Mesh, bool Const>
EdgeHandle<Mesh, Const> HalfedgeHandle<Mesh, Const>::edge() const {
  return EdgeHandle<Mesh, Const>{mesh->he_edge(id), mesh};
}

template <class Mesh, bool Const>
HalfedgeHandle<Mesh, Const> HalfedgeHandle<Mesh, Const>::sibling() const {
  return HalfedgeHandle<Mesh, Const>{mesh->he_sibling(id), mesh};
}

template <class Mesh, bool Const>
HalfedgeHandle<Mesh, Const> HalfedgeHandle<Mesh, Const>::incoming_next() const {
  return HalfedgeHandle<Mesh, Const>{mesh->he_incoming_next(id), mesh};
}

template <class Mesh, bool Const>
HalfedgeHandle<Mesh, Const> HalfedgeHandle<Mesh, Const>::twin() const {
  const VertexId to_vid = data().vertex;
  HalfedgeHandle<Mesh, Const> sib = sibling();
  while (sib.data().vertex == to_vid) {
    sib = sib.sibling();
  }
  return sib;
}

template <class Mesh, bool Const>
HalfedgeHandle<Mesh, Const> EdgeHandle<Mesh, Const>::halfedge() const {
  return HalfedgeHandle<Mesh, Const>{mesh->e_halfedge(id), mesh};
}

template <class Mesh, bool Const>
auto EdgeHandle<Mesh, Const>::halfedges() const -> HalfedgeCycleRange<Mesh, Const, &Mesh::he_sibling> {
  return HalfedgeCycleRange<Mesh, Const, &Mesh::he_sibling>{mesh, mesh->e_halfedge(id)};
}

template <class Mesh, bool Const>
std::array<VertexHandle<Mesh, Const>, 2> EdgeHandle<Mesh, Const>::vertices() const {
  const HalfedgeHandle<Mesh, Const> he = halfedge();
  const HalfedgeHandle<Mesh, Const> prev = he.prev();
  return {prev.to(), he.to()};
}

}  // namespace gpf
