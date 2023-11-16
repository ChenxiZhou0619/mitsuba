// clang-format off
#pragma once
#include <mitsuba/mitsuba.h>
#include <cmath>

// clang-format on

MTS_NAMESPACE_BEGIN

#define MAX_SH_ORDER 64

class SphericalHarmonic {
public:
  SphericalHarmonic() = delete;

  SphericalHarmonic(uint32_t maxOrder) {
    if (maxOrder > MAX_SH_ORDER) {
      printf("Exceed maximum sh order\n");
      exit(1);
    }
    m_maxOrder            = maxOrder;
    uint32_t nCofficients = (maxOrder + 1) * (maxOrder + 1);
    m_cofficients         = std::vector<Float>(nCofficients);
  };

  static inline Float Basis(uint32_t l, int m, Float theta, Float phi) {
    Float legendre = std::sph_legendre(l, std::abs(m), theta);
    Float phiTerm  = (m >= 0) ? std::cos(m * phi) : (std::sin(-m * phi));
    phiTerm *= (m == 0) ? 1 : std::sqrt(2);
    return phiTerm * legendre;
  }

  void addSample(Vector v, Float f);

  Float eval(Vector v) const;

private:
  uint32_t           m_maxOrder;
  std::vector<Float> m_cofficients;
};

inline void SphericalHarmonic::addSample(Vector v, Float f) {
  int   offset = 0;
  Float theta  = math::safe_acos(v.y);
  Float phi    = std::atan2(v.x, -v.z);
  if (phi < 0) phi += 2 * M_PI;

  for (int l = 0; l <= m_maxOrder; ++l) {
    for (int m = -l; m <= l; ++m) {
      Float a = Basis(l, m, theta, phi);
      m_cofficients[offset++] += a * f;
    }
  }
}

inline Float SphericalHarmonic::eval(Vector v) const {
  int   offset = 0;
  Float val    = .0f;

  Float theta = math::safe_acos(v.y);
  Float phi   = std::atan2(v.x, -v.z);
  if (phi < 0) phi += 2 * M_PI;

  for (int l = 0; l <= m_maxOrder; ++l) {
    for (int m = -l; m <= l; ++m) {
      val += m_cofficients[offset++] * Basis(l, m, theta, phi);
    }
  }

  return val;
}

MTS_NAMESPACE_END
