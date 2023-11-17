#include "spectralintegrator.h"
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/sensor.h>
MTS_NAMESPACE_BEGIN

Spectrum SpectralMonteCarloIntegrator::Li(const RayDifferential &r,
                                          RadianceQueryRecord   &rRec) const {
  //
  Log(EError, "SpectralMonteCarloIntegretor shouldn'd call this function");
  exit(1);
}

void SpectralMonteCarloIntegrator::renderBlock(
    const Scene *scene, const Sensor *sensor, Sampler *sampler,
    ImageBlock *block, const bool &stop,
    const std::vector<TPoint2<uint8_t>> &points) const {

  Float diffScaleFactor = 1.0f / std::sqrt((Float)sampler->getSampleCount());

  bool needsApertureSample = sensor->needsApertureSample();
  bool needsTimeSample     = sensor->needsTimeSample();
  bool needsSpectrumSample = true;

  RadianceQueryRecord rRec(scene, sampler);
  Point2              apertureSample(0.5f);
  Float               timeSample   = 0.5f;
  Float               lambdaSample = .5f;
  RayDifferential     sensorRay;

  block->clear();

  uint32_t queryType = RadianceQueryRecord::ESensorRay;

  if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we
                                         don't have to */
    queryType &= ~RadianceQueryRecord::EOpacity;

  for (size_t i = 0; i < points.size(); ++i) {
    Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
    if (stop) break;

    sampler->generate(offset);

    for (size_t j = 0; j < sampler->getSampleCount(); j++) {
      rRec.newQuery(queryType, sensor->getMedium());
      Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

      if (needsApertureSample) apertureSample = rRec.nextSample2D();
      if (needsTimeSample) timeSample = rRec.nextSample1D();
      if (needsSpectrumSample) lambdaSample = rRec.nextSample1D();

      Spectrum sensorWeight = sensor->sampleRayDifferential(
          sensorRay, samplePos, apertureSample, timeSample);
      spectral::SampledWavelengths lambdas =
          spectral::SampledWavelengths::SampleUniform(
              lambdaSample); // TODO This could be configured

      sensorRay.scaleDifferential(diffScaleFactor);
      spectral::SampledSpectrum spec = Li(sensorRay, rRec, lambdas);

      //* Turn SampledSpectrum to RGB
      Float rgb[3];
      spectral::SampledSpectrumToRGB(spec, lambdas, rgb);

      Spectrum rgbLi = sensorWeight * Spectrum{rgb};

      block->put(samplePos, rgbLi, rRec.alpha);
      sampler->advance();
    }
  }
}

MTS_IMPLEMENT_CLASS(MonteCarloIntegrator, true, MonteCarloIntegrator)

MTS_NAMESPACE_END
