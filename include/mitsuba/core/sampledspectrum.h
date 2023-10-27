#pragma once

#include <array>
#include <mitsuba/mitsuba.h>
MTS_NAMESPACE_BEGIN

#define NSpectrumSamples 4

//! SampledSpectrum doesn't store the lambda information
//! Make sure the different SampledSpectrum objects have same sampled wavelengths

class Spectrum;

class SampledWavelengths {
public:
  static SampledWavelengths SampleUniform(Float u);

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

  Spectrum toRGB(const SampledWavelengths &lambda) const;

private:
  std::array<Float, NSpectrumSamples> values;
};

SampledSpectrum operator*(Float f, const SampledSpectrum &spec);

SampledSpectrum RGBIlluminantToSampledSpectrum(const Spectrum           &rgb,
                                               const SampledWavelengths &lambda);
SampledSpectrum RGBAlbedoToSampledSpectrum(const Spectrum &rgb, const SampledWavelengths &lambda);

MTS_NAMESPACE_END