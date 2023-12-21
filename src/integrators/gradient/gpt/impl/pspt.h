#pragma once

#include "../gradient_pt.h"

MTS_NAMESPACE_BEGIN

/// Gradient Path Tracer implementation with Primary Space shift mapping
class PSGradientPathTracer : public IGradientPathTracer {
public:
  PSGradientPathTracer(const Scene *scene, const Sensor *sensor,
                       Sampler *sampler, GPTWorkResult *block,
                       const GradientPathTracerConfig *config)
      : IGradientPathTracer(scene, sensor, sampler, block, config) {}

  /// Evaluates a sample at the given position.
  ///
  /// Outputs direct radiance to be added on top of the final image, the
  /// throughput to the central pixel, gradients to all neighbors, and
  /// throughput contribution to the neighboring pixels.
  virtual void evaluatePoint(RadianceQueryRecord &rRec,
                             const Point2        &samplePosition,
                             const Point2 &apertureSample, Float timeSample,
                             Float     differentialScaleFactor,
                             Spectrum &out_very_direct,
                             Spectrum &out_throughput, Spectrum *out_gradients,
                             Spectrum *out_neighborThroughputs) override;

  /// Constructs a sequence of base paths and shifts them into offset paths,
  /// evaluating their throughputs and differences.
  ///
  /// This is the core of the rendering algorithm.
  virtual void evaluate(RayState &main, RayState *shiftedRays,
                        int secondaryCount, Spectrum &out_veryDirect) override;
};

MTS_NAMESPACE_END