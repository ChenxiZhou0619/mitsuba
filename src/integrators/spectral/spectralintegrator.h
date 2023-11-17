// clang-format off
#pragma once
#include <mitsuba/render/integrator.h>
#include "spectrum/spectrum.h"
// clang-format on

MTS_NAMESPACE_BEGIN

class SpectralMonteCarloIntegrator : public MonteCarloIntegrator {
public:
  void serialize(Stream *stream, InstanceManager *manager) const override {
    MonteCarloIntegrator::serialize(stream, manager);
    Log(EError, "SpectralMonteCarloIntegrator serialization not support");
  }

  virtual void
  renderBlock(const Scene *scene, const Sensor *sensor, Sampler *sampler,
              ImageBlock *block, const bool &stop,
              const std::vector<TPoint2<uint8_t>> &points) const override;

  //! This function is deprecated since no lambda information is provided during
  //! rendering
  virtual Spectrum Li(const RayDifferential &r,
                      RadianceQueryRecord   &rRec) const override final;

  //! This should be called
  virtual spectral::SampledSpectrum
  Li(const RayDifferential &r, RadianceQueryRecord &rRec,
     spectral::SampledWavelengths &lambdas) const = 0;

  MTS_DECLARE_CLASS()

protected:
  SpectralMonteCarloIntegrator(const Properties &props)
      : MonteCarloIntegrator(props) {
    //
  }

  SpectralMonteCarloIntegrator(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {
    Log(EError, "SpectralMonteCarloIntegrator serialization not support");
  }

  virtual ~SpectralMonteCarloIntegrator() {}
};

MTS_NAMESPACE_END