// clang-format off
#pragma once
#include <mitsuba/mitsuba.h>
// clang-format on

MTS_NAMESPACE_BEGIN

class SphericalHarmonic {
public:
  SphericalHarmonic() = delete;

  SphericalHarmonic(uint32_t maxOrder) : m_maxOrder(maxOrder) {
    m_cofficients = std::vector<Float>((maxOrder + 1) * (maxOrder + 1));
  };

  //

private:
  uint32_t           m_maxOrder;
  std::vector<Float> m_cofficients;
};
MTS_NAMESPACE_END
