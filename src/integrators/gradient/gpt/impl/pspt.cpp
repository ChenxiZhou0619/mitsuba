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
  mainRay.rRec = rRec;
  mainRay.rRec.its = rRec.its;

  /// Offset path
  RayState shiftedRays[4];
  static const Vector2 pixelShifts[4] = {
      Vector2(1.0f, 0.0f), Vector2(0.0f, 1.0f), Vector2(-1.0f, 0.0f),
      Vector2(0.0f, -1.0f)};

  for (int i = 0; i < 4; ++i) {
    shiftedRays[i].throughput = m_sensor->sampleRayDifferential(
        shiftedRays[i].ray, samplePosition + pixelShifts[i], apertureSample,
        timeSample);
    shiftedRays[i].ray.scaleDifferential(differentialScaleFactor);
    shiftedRays[i].rRec = rRec;
    shiftedRays[i].rRec.its = rRec.its;
  }

  // Evaluate the gradients. The actual algorithm happens here.
  Spectrum very_direct = Spectrum(0.0f);
  evaluate(mainRay, shiftedRays, 4, very_direct);

  // Output results.
  out_very_direct = very_direct;
  out_throughput = mainRay.radiance;

  for (int i = 0; i < 4; i++) {
    out_gradients[i] = shiftedRays[i].gradient;
    out_neighborThroughputs[i] = shiftedRays[i].radiance;
  }
}

void PSGradientPathTracer::evaluate(RayState &main, RayState *shiftedRays,
                                    int secondaryCount,
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
    // First hit is not in the scene so can't continue. Also there there are
    // no paths to shift.

    // Add potential very direct light from the environment as gradients are
    // not used for that.
    if (main.rRec.type & RadianceQueryRecord::EEmittedRadiance) {
      out_veryDirect += main.throughput * scene->evalEnvironment(main.ray);
    }
    return;
  }

  // Add very direct light from non-environment.
  {
    // Include emitted radiance if requested.
    if (main.rRec.its.isEmitter() &&
        (main.rRec.type & RadianceQueryRecord::EEmittedRadiance)) {
      out_veryDirect += main.throughput * main.rRec.its.Le(-main.ray.d);
    }
  }

  // If no intersection of an offset ray could be found, its offset paths can
  // not be generated.
  for (int i = 0; i < secondaryCount; ++i) {
    RayState &shifted = shiftedRays[i];
    if (!shifted.rRec.its.isValid()) {
      shifted.alive = false;
    }
  }

  // Strict normals check to produce the same results as bidirectional methods
  // when normal mapping is used.
  if (m_config->m_strictNormals) {
    // If 'strictNormals'=true, when the geometric and shading normals
    // classify the incident direction to the same side, then the main path is
    // still good.
    if (dot(main.ray.d, main.rRec.its.geoFrame.n) *
            Frame::cosTheta(main.rRec.its.wi) >=
        0) {
      // This is an impossible base path.
      return;
    }

    for (int i = 0; i < secondaryCount; ++i) {
      RayState &shifted = shiftedRays[i];

      if (dot(shifted.ray.d, shifted.rRec.its.geoFrame.n) *
              Frame::cosTheta(shifted.rRec.its.wi) >=
          0) {
        // This is an impossible offset path.
        shifted.alive = false;
      }
    }
  }

  // Main path tracing loop.
  main.rRec.depth = 1;

  while (main.rRec.depth < m_config->m_maxDepth || m_config->m_maxDepth < 0) {
    // Strict normals check to produce the same results as bidirectional
    // methods when normal mapping is used. If 'strictNormals'=true, when the
    // geometric and shading normals classify the incident direction to the
    // same side, then the main path is still good.
    if (m_config->m_strictNormals) {
      if (dot(main.ray.d, main.rRec.its.geoFrame.n) *
              Frame::cosTheta(main.rRec.its.wi) >=
          0) {
        // This is an impossible main path, and there are no more paths to
        // shift.
        return;
      }

      for (int i = 0; i < secondaryCount; ++i) {
        RayState &shifted = shiftedRays[i];

        if (dot(shifted.ray.d, shifted.rRec.its.geoFrame.n) *
                Frame::cosTheta(shifted.rRec.its.wi) >=
            0) {
          // This is an impossible offset path.
          shifted.alive = false;
        }
      }
    }

    //------ Next Event Estimation ------
    {
      // 1. estimate the main ray

      const BSDF *main_bsdf = main.rRec.its.getBSDF(main.ray);

      if (main.rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
          main_bsdf->getType() &
              BSDF::ESmooth /* Apply Nee when bsdf is not delta*/) {
        DirectSamplingRecord dRec{main.rRec.its};
        const Point2 light_sample = main.rRec.nextSample2D();

        auto [Li, visiable] =
            m_scene->sampleEmitterDirectVisible(dRec, light_sample);

        const Emitter *emitter = static_cast<const Emitter *>(dRec.object);
        SAssert(emitter != nullptr);

        BSDFSamplingRecord main_bsdf_bRec{
            main.rRec.its, main.rRec.its.toLocal(dRec.d), ERadiance};
        Spectrum main_bsdf_value = main_bsdf->eval(main_bsdf_bRec);

        Float pdf_bsdf =
            (emitter->isOnSurface() && dRec.measure == ESolidAngle && visiable)
                ? main_bsdf->pdf(main_bsdf_bRec)
                : .0f;

        Float misw = dRec.pdf / (dRec.pdf + pdf_bsdf);

        // Li =  throughput * misw * bsdf_value * Li / pdf
        Spectrum main_radiance = main.throughput * misw * main_bsdf_value * Li;

        main.addRadiance(main_radiance, 1.f);

        // 2. estimate the shifted rays
        for (int i = 0; i < 4; ++i) {
          RayState &shifted = shiftedRays[i];

          bool shift_successful = shifted.alive;
          if (shift_successful) {
            // apply primary space shift mapping
            // i.e. reuse light_sample

            DirectSamplingRecord shifted_dRec{shifted.rRec.its};
            auto [shifted_Li, shifted_visiable] =
                m_scene->sampleEmitterDirectVisible(shifted_dRec, light_sample);

            const Emitter *emitter =
                static_cast<const Emitter *>(shifted_dRec.object);
            SAssert(emitter != nullptr);

            const BSDF *shifted_bsdf = shifted.rRec.its.getBSDF();
            BSDFSamplingRecord shifted_bsdf_bRec{
                shifted.rRec.its, shifted.rRec.its.toLocal(shifted_dRec.d),
                ERadiance};
            Spectrum shifted_bsdf_value = shifted_bsdf->eval(shifted_bsdf_bRec);

            Float shifted_pdf_bsdf =
                (emitter->isOnSurface() && dRec.measure == ESolidAngle &&
                 shifted_visiable)
                    ? shifted_bsdf->pdf(shifted_bsdf_bRec)
                    : .0f;
            Float shifted_misw =
                shifted_dRec.pdf / (shifted_dRec.pdf + shifted_pdf_bsdf);

            // shifted_Li = throughput * misw * jacobian * bsdf_value * Li / pdf
            Spectrum shifted_radiance = shifted.throughput * shifted_misw *
                                        shifted_bsdf_value * shifted_Li;

            shifted.addRadiance(shifted_radiance, 1.f);
          }
        }
      }
    }

    //------ BSDF Sampling Estimation -----
    {
      Point2 bsdf_sample = main.rRec.nextSample2D();
      const BSDF *main_bsdf = main.rRec.its.getBSDF(main.ray);

      Float main_bsdf_pdf;
      BSDFSamplingRecord main_bRec{main.rRec.its, main.rRec.sampler, ERadiance};
      Spectrum main_bsdf_weight =
          main_bsdf->sample(main_bRec, main_bsdf_pdf, bsdf_sample);

      if (main_bsdf_weight.isZero())
        break;

      DirectSamplingRecord main_dRec{main.rRec.its};
      bool main_hitter_emitter = false;

      const Vector main_wo = main.rRec.its.toWorld(main_bRec.wo);
      Intersection previous_main_its = main.rRec.its;
      main.ray = Ray{main.rRec.its.p, main_wo, main.ray.time};

      Spectrum Le{.0f};

      if (scene->rayIntersect(main.ray, main.rRec.its)) {
        if (main.rRec.its.isEmitter() /* If hit the emitter*/) {
          Le = main.rRec.its.Le(-main.ray.d);
          main_dRec.setQuery(main.ray, main.rRec.its);
          main_hitter_emitter = true;
        }
      } else {
        const Emitter *env = scene->getEnvironmentEmitter();

        if (env) {
          Le = env->evalEnvironment(main.ray);
          if (!env->fillDirectSamplingRecord(main_dRec, main.ray))
            break;
          main_hitter_emitter = true;
        }
      }

      main.throughput *= main_bsdf_weight;
      if (main_hitter_emitter &&
          (main.rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
        const Float main_lum_pdf = (!(main_bRec.sampledType & BSDF::EDelta))
                                       ? scene->pdfEmitterDirect(main_dRec)
                                       : 0;
        Float misw = main_bsdf_pdf / (main_bsdf_pdf + main_lum_pdf);
        Spectrum main_radiance = main.throughput * misw * Le;
        main.addRadiance(main_radiance, 1.f);
      }
    }
  }
}

MTS_NAMESPACE_END