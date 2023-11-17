#pragma once
#include <array>
#include <mitsuba/mitsuba.h>
using Float = mitsuba::Float;

namespace spectral {

constexpr static size_t NSpectrumSamples = 4;
constexpr static Float  LambdaMin        = 360;
constexpr static Float  LambdaMax        = 830;

//! SampledSpectrum doesn't store the lambda information
//! Make sure the different SampledSpectrum objects have same sampled
//! wavelengths
class SampledWavelengths {
public:
  static SampledWavelengths SampleUniform(Float u);

  Float operator[](int entry) const;

private:
  std::array<Float, NSpectrumSamples> wavelengths, pdfs;
};

class SampledSpectrum {
public:
  explicit SampledSpectrum(Float v = 0.0f);

  Float operator[](int entry) const;

  Float &operator[](int entry);

  SampledSpectrum operator+(const SampledSpectrum &spec) const;

  SampledSpectrum &operator+=(const SampledSpectrum &spec);

  SampledSpectrum operator-(const SampledSpectrum &spec) const;

  SampledSpectrum &operator-=(const SampledSpectrum &spec);

  SampledSpectrum operator*(const SampledSpectrum &spec) const;

  SampledSpectrum &operator*=(const SampledSpectrum &spec);

  SampledSpectrum operator*(Float f) const;

  SampledSpectrum &operator*=(Float f);

  SampledSpectrum operator/(const SampledSpectrum &spec) const;

  SampledSpectrum &operator/=(const SampledSpectrum &spec);

  SampledSpectrum operator/(Float f) const;

  SampledSpectrum &operator/=(Float f);

  bool operator==(const SampledSpectrum &spec) const;

  bool operator!=(const SampledSpectrum &spec) const;

  explicit operator bool() const;

  bool isNaN() const;

  bool isValid() const;

  Float max() const;

  Float min() const;

  Float average() const;

private:
  std::array<Float, NSpectrumSamples> values;
};

/**
 * RGBAlbedoSpectrum, represented by a sigmoid polynomial
 */

class RGBAlbedoSpectrum {
public:
  explicit RGBAlbedoSpectrum(Float v = .0f);

  explicit RGBAlbedoSpectrum(Float r, Float g, Float b);

  explicit RGBAlbedoSpectrum(const mitsuba::Spectrum &rgb);

  inline Float operator()(Float lambda) const;

  SampledSpectrum sample(const SampledWavelengths &lambdas) const;

private:
  Float m_coeffs[3];
};

class RGBToSpectrumTable {
public:
  static constexpr int res = 64;
  using CoefficientArray   = Float[3][res][res][res][3];

  RGBToSpectrumTable(const Float *zNodes, const CoefficientArray *coeffs)
      : m_zNodes(zNodes), m_coeffs(coeffs) {}

  static void Init();

  void RGBToPolyCoeff(const Float *rgb, Float *coeffs) const;

public:
  // support color space
  static std::unique_ptr<RGBToSpectrumTable> sRGB;

private:
  const Float            *m_zNodes;
  const CoefficientArray *m_coeffs;
};

// TODO consider the color space
void SampledSpectrumToRGB(const SampledSpectrum    &spec,
                          const SampledWavelengths &lambdas, Float *rgb);

SampledSpectrum operator*(Float f, const SampledSpectrum &spec);

enum RGBType {
  EIlluminantRGB  = 1 << 0,
  EReflectanceRGB = 1 << 1,
  EUnboundedRGB   = 1 << 2
};

// TODO consider the color space
SampledSpectrum RGBToSampledSpectrum(const mitsuba::Spectrum  &rgb,
                                     const SampledWavelengths &lambdas,
                                     RGBType                   type);

namespace math {
template <typename Predicate>
inline size_t FindInterval(size_t sz, const Predicate &pred) {
  using ssize_t = std::make_signed_t<size_t>;
  ssize_t size = (ssize_t)sz - 2, first = 1;
  while (size > 0) {
    // Evaluate predicate at midpoint and update _first_ and _size_
    size_t half = (size_t)size >> 1, middle = first + half;
    bool   predResult = pred(middle);
    first             = predResult ? middle + 1 : first;
    size              = predResult ? size - (half + 1) : half;
  }
  return (size_t)std::clamp<ssize_t>((ssize_t)first - 1, 0, sz - 2);
}

inline Float Lerp(Float x, Float a, Float b) { return (1 - x) * a + x * b; }

}; // namespace math

} // namespace spectral
