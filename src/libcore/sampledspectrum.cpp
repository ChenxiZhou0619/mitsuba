#include <mitsuba/core/sampledspectrum.h>

MTS_NAMESPACE_BEGIN

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

//! Convert point-sampled spectra to RGB coefficients
Spectrum SampledSpectrum::toRGB(const SampledWavelengths &lambda) const {
  //
}

MTS_NAMESPACE_END