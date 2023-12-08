#include "pathgrid.h"

MTS_NAMESPACE_BEGIN

namespace pathgrid {

PathGrid::PathGrid(uint32_t max_size)
    : m_max_size(max_size), m_vertices(max_size) {

  m_wis.reserve(max_size);
  m_Lis.reserve(max_size);
}

} // namespace pathgrid

MTS_NAMESPACE_END