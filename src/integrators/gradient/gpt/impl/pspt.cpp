#include "pspt.h"
#include "shiftmap.h"

MTS_NAMESPACE_BEGIN

void PSGradientPathTracer::evaluatePoint(
    RadianceQueryRecord &rRec, const Point2 &samplePosition,
    const Point2 &apertureSample, Float timeSample,
    Float differentialScaleFactor, Spectrum &out_very_direct,
    Spectrum &out_throughput, Spectrum *out_gradients,
    Spectrum *out_neighborThroughputs) {

  /// Base path
  RayState mainRay;
  mainRay.throughput = m_sensor->sampleRayDifferential(
      mainRay.ray, samplePosition, apertureSample, timeSample);
  mainRay.ray.scaleDifferential(differentialScaleFactor);
  mainRay.rRec     = rRec;
  mainRay.rRec.its = rRec.its;

  /// Offset path
  RayState             shiftedRays[4];
  static const Vector2 pixelShifts[4] = {
      Vector2(1.0f, 0.0f), Vector2(0.0f, 1.0f), Vector2(-1.0f, 0.0f),
      Vector2(0.0f, -1.0f)};

  for (int i = 0; i < 4; ++i) {
    shiftedRays[i].throughput = m_sensor->sampleRayDifferential(
        shiftedRays[i].ray, samplePosition + pixelShifts[i], apertureSample,
        timeSample);
    shiftedRays[i].ray.scaleDifferential(differentialScaleFactor);
    shiftedRays[i].rRec     = rRec;
    shiftedRays[i].rRec.its = rRec.its;
  }

  // Evaluate the gradients. The actual algorithm happens here.
  Spectrum very_direct = Spectrum(0.0f);
  evaluate(mainRay, shiftedRays, 4, very_direct);

  // Output results.
  out_very_direct = very_direct;
  out_throughput  = mainRay.radiance;

  for (int i = 0; i < 4; i++) {
    out_gradients[i]           = shiftedRays[i].gradient;
    out_neighborThroughputs[i] = shiftedRays[i].radiance;
  }
}

void PSGradientPathTracer::evaluate(RayState &main, RayState *shiftedRays,
                                    int       secondaryCount,
                                    Spectrum &out_veryDirect) {
  const Scene *scene = main.rRec.scene;

  // Base path intersection
  main.rRec.rayIntersect(main.ray);
  main.ray.mint = Epsilon;

  for (int i = 0; i < secondaryCount; ++i) {
    RayState &shifted = shiftedRays[i];
    shifted.rRec.rayIntersect(shifted.ray);
    shifted.ray.mint = Epsilon;
  }

  if (!main.rRec.its.isValid()) {
    //
  }
}

MTS_NAMESPACE_END