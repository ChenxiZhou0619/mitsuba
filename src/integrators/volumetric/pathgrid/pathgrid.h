// clang-format off
#pragma once
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>
#include "path.h"
#include "nanoflann.hpp"
// clang-format on

MTS_NAMESPACE_BEGIN

namespace pathgrid {

struct VertexInfo {
  Point    p;
  Vector   wi;
  Spectrum contrib;
};

class PathGrid {
public:
  using kdtree = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<Float, PathGrid>, PathGrid, 3>;

  PathGrid(int max_size) : m_max_size(max_size) {
    m_pts.reserve(max_size);
    m_wis.reserve(max_size);
    m_contrib.reserve(max_size);
  }

  virtual ~PathGrid() { delete m_kdtree; }

  void addStorage(const PathStorage &storage) {
    int i = 0;
    while (m_pts.size() < m_max_size) {
      m_pts.emplace_back(storage.ps[i]);
      m_wis.emplace_back(storage.wis[i]);
      m_contrib.emplace_back(storage.contribs[i]);
      ++i;
    }
  }

  void init() { m_kdtree = new kdtree(3, *this, {10}); }

  inline size_t kdtree_get_point_count() const { return m_pts.size(); }

  inline Float kdtree_get_pt(const size_t idx, const size_t dim) const {
    return m_pts[idx][dim];
  }

  template <class BBox> bool kdtree_get_bbox(BBox & /* bb */) const {
    return false;
  }

  int searchKNN(const Float *query_pt, int num_result, uint32_t *ret_idx,
                Float *ret_distance_sqr) const {
    num_result =
        m_kdtree->knnSearch(query_pt, num_result, ret_idx, ret_distance_sqr);
    return num_result;
  }

  std::tuple<Point, Vector, Spectrum> getData(int idx) const {
    return {m_pts[idx], m_wis[idx], m_contrib[idx]};
  }

protected:
  int                   m_max_size;
  std::vector<Point>    m_pts;
  std::vector<Vector>   m_wis;
  std::vector<Spectrum> m_contrib;

  kdtree *m_kdtree = nullptr;
};

} // namespace pathgrid

MTS_NAMESPACE_END
