// clang-format off
#pragma once
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>
// clang-format on

MTS_NAMESPACE_BEGIN

namespace pathgrid {

struct VertexSet {

  struct Position {
    Float x, y, z;
  };

  std::vector<Position> pts;

public:
  VertexSet(uint32_t max_size);
};

class PathGrid {
public:
  PathGrid(uint32_t max_size);

  void searchKNN(uint32_t idx[], uint32_t K) const;

protected:
  uint32_t              m_max_size;
  VertexSet             m_vertices;
  std::vector<Vector>   m_wis;
  std::vector<Spectrum> m_Lis;
};

} // namespace pathgrid

MTS_NAMESPACE_END
