#include "spectrum.h"

namespace spectral {

// TODO
SampledWavelengths SampledWavelengths::SampleUniform(Float u) {
  //
}

Float SampledWavelengths::operator[](int entry) const {
  return wavelengths[entry];
}

SampledSpectrum::SampledSpectrum(Float v) {
  for (int i = 0; i < NSpectrumSamples; ++i)
    values[i] = v;
}

Float SampledSpectrum::operator[](int entry) const { return values[entry]; }

Float &SampledSpectrum::operator[](int entry) { return values[entry]; }

SampledSpectrum SampledSpectrum::operator+(const SampledSpectrum &spec) const {
  SampledSpectrum s;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    s[i] = values[i] + spec[i];
  }
  return s;
}

SampledSpectrum &SampledSpectrum::operator+=(const SampledSpectrum &spec) {
  for (int i = 0; i < NSpectrumSamples; ++i)
    values[i] += spec[i];
  return *this;
}

SampledSpectrum SampledSpectrum::operator-(const SampledSpectrum &spec) const {
  SampledSpectrum s;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    s[i] = values[i] - spec[i];
  }
  return s;
}

SampledSpectrum &SampledSpectrum::operator-=(const SampledSpectrum &spec) {
  for (int i = 0; i < NSpectrumSamples; ++i)
    values[i] -= spec[i];
  return *this;
}

SampledSpectrum SampledSpectrum::operator*(const SampledSpectrum &spec) const {
  SampledSpectrum s;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    s[i] = values[i] * spec[i];
  }
  return s;
}

SampledSpectrum SampledSpectrum::operator*(Float f) const {
  SampledSpectrum s;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    s[i] = values[i] * f;
  }
  return s;
}

SampledSpectrum &SampledSpectrum::operator*=(const SampledSpectrum &spec) {
  for (int i = 0; i < NSpectrumSamples; ++i)
    values[i] *= spec[i];
  return *this;
}

SampledSpectrum &SampledSpectrum::operator*=(Float f) {
  for (int i = 0; i < NSpectrumSamples; ++i)
    values[i] *= f;
  return *this;
}

SampledSpectrum SampledSpectrum::operator/(const SampledSpectrum &spec) const {
  SampledSpectrum s;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    s[i] = values[i] / spec[i];
  }
  return s;
}

SampledSpectrum SampledSpectrum::operator/(Float f) const {
  SampledSpectrum s;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    s[i] = values[i] / f;
  }
  return s;
}

SampledSpectrum &SampledSpectrum::operator/=(const SampledSpectrum &spec) {
  for (int i = 0; i < NSpectrumSamples; ++i)
    values[i] /= spec[i];
  return *this;
}

SampledSpectrum &SampledSpectrum::operator/=(Float f) {
  for (int i = 0; i < NSpectrumSamples; ++i)
    values[i] /= f;
  return *this;
}

bool SampledSpectrum::operator==(const SampledSpectrum &spec) const {
  for (int i = 0; i < NSpectrumSamples; ++i)
    if (values[i] != spec[i]) return false;
  return true;
}

bool SampledSpectrum::operator!=(const SampledSpectrum &spec) const {
  for (int i = 0; i < NSpectrumSamples; ++i)
    if (values[i] != spec[i]) return true;
  return false;
}

SampledSpectrum::operator bool() const {
  for (int i = 0; i < NSpectrumSamples; ++i)
    if (values[i] != .0f) return true;
  return false;
}

bool SampledSpectrum::isNaN() const {
  for (int i = 0; i < NSpectrumSamples; ++i)
    if (std::isnan(values[i])) return true;
  return false;
}

bool SampledSpectrum::isValid() const {
  for (int i = 0; i < NSpectrumSamples; ++i)
    if (!std::isfinite(values[i])) return false;
  return true;
}

Float SampledSpectrum::max() const {
  Float m = values[0];
  for (int i = 1; i < NSpectrumSamples; ++i)
    m = std::max(m, values[i]);
  return m;
}

Float SampledSpectrum::min() const {
  Float m = values[0];
  for (int i = 1; i < NSpectrumSamples; ++i)
    m = std::min(m, values[i]);
  return m;
}

Float SampledSpectrum::average() const {
  Float mean = .0f;
  for (int i = 0; i < NSpectrumSamples; ++i)
    mean += values[i];
  return mean / NSpectrumSamples;
}

SampledSpectrum operator*(Float f, const SampledSpectrum &spec) {
  SampledSpectrum s;
  for (int i = 0; i < NSpectrumSamples; ++i)
    s[i] = f * spec[i];
  return s;
}

extern const int                                  sRGBToSpectrumTable_Res;
extern const Float                                sRGBToSpectrumTable_Scale[64];
extern const RGBToSpectrumTable::CoefficientArray sRGBToSpectrumTable_Data;

void RGBToSpectrumTable::Init() {
  sRGB = std::make_unique<RGBToSpectrumTable>(sRGBToSpectrumTable_Scale,
                                              &sRGBToSpectrumTable_Data);
}

void RGBToSpectrumTable::RGBToPolyCoeff(const Float *rgb, Float *coeffs) const {

  // handle uniform rgb
  if (rgb[0] == rgb[1] && rgb[1] == rgb[2]) {
    coeffs[0] = .0f;
    coeffs[1] = .0f;
    coeffs[2] = (rgb[0] - .5f) / std::sqrt(rgb[0] * (1 - rgb[0]));
    return;
  }

  int   maxc = (rgb[0] > rgb[1]) ? ((rgb[0] > rgb[2]) ? 0 : 2)
                                 : ((rgb[1] > rgb[2]) ? 1 : 2);
  float z    = rgb[maxc];
  float x    = rgb[(maxc + 1) % 3] * (res - 1) / z;
  float y    = rgb[(maxc + 2) % 3] * (res - 1) / z;

  // Compute integer indices and offsets for coefficient interpolation
  int xi = std::min((int)x, res - 2), yi = std::min((int)y, res - 2),
      zi   = spectral::math::FindInterval(res,
                                          [&](int i) { return m_zNodes[i] < z; });
  Float dx = x - xi, dy = y - yi,
        dz = (z - m_zNodes[zi]) / (m_zNodes[zi + 1] - m_zNodes[zi]);

  for (int i = 0; i < 3; ++i) {
    // Define _co_ lambda for looking up sigmoid polynomial coefficients
    auto co = [&](int dx, int dy, int dz) {
      return (*m_coeffs)[maxc][zi + dz][yi + dy][xi + dx][i];
    };

    coeffs[i] = spectral::math::Lerp(
        dz,
        spectral::math::Lerp(
            dy, spectral::math::Lerp(dx, co(0, 0, 0), co(1, 0, 0)),
            spectral::math::Lerp(dx, co(0, 1, 0), co(1, 1, 0))),
        spectral::math::Lerp(
            dy, spectral::math::Lerp(dx, co(0, 0, 1), co(1, 0, 1)),
            spectral::math::Lerp(dx, co(0, 1, 1), co(1, 1, 1))));
  }
}

RGBAlbedoSpectrum::RGBAlbedoSpectrum(Float v) {
  Float rgb[3] = {v, v, v};
  RGBToSpectrumTable::sRGB->RGBToPolyCoeff(rgb, m_coeffs);
}
RGBAlbedoSpectrum::RGBAlbedoSpectrum(Float r, Float g, Float b) {
  Float rgb[3] = {r, g, b};
  RGBToSpectrumTable::sRGB->RGBToPolyCoeff(rgb, m_coeffs);
}

RGBAlbedoSpectrum::RGBAlbedoSpectrum(const mitsuba::Spectrum &spec) {
  Float rgb[3] = {spec[0], spec[1], spec[2]};
  RGBToSpectrumTable::sRGB->RGBToPolyCoeff(rgb, m_coeffs);
}

Float RGBAlbedoSpectrum::operator()(Float lambda) const {
  auto SigmoidPoly = [&](Float lambda) {
    Float poly =
        m_coeffs[0] * lambda * lambda + m_coeffs[1] * lambda + m_coeffs[2];
    if (std::isinf(poly)) {
      return (poly > 0) ? 1.f : 0.f;
    }
    return .5f + poly / (2.f * std::sqrt(1 + poly * poly));
  };

  return SigmoidPoly(lambda);
}

SampledSpectrum
RGBAlbedoSpectrum::sample(const SampledWavelengths &lambdas) const {
  SampledSpectrum spec;

  spec[0] = (*this)(lambdas[0]);
  spec[1] = (*this)(lambdas[1]);
  spec[2] = (*this)(lambdas[2]);
  spec[3] = (*this)(lambdas[3]);

  return spec;
}

}; // namespace spectral
