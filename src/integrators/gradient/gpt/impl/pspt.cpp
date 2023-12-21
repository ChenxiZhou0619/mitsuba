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
    // First hit is not in the scene so can't continue. Also there there are
    // no paths to shift.

    // Add potential very direct light from the environment as gradients are
    // not used for that.
    if (main.rRec.type & RadianceQueryRecord::EEmittedRadiance) {
      out_veryDirect += main.throughput * scene->evalEnvironment(main.ray);
    }

    // SLog(EInfo, "Main ray(%d): First hit not in scene.", rayCount);
    return;
  }

  // Add very direct light from non-environment.
  {
    // Include emitted radiance if requested.
    if (main.rRec.its.isEmitter() &&
        (main.rRec.type & RadianceQueryRecord::EEmittedRadiance)) {
      out_veryDirect += main.throughput * main.rRec.its.Le(-main.ray.d);
    }

    // Include radiance from a subsurface scattering model if requested. Note:
    // Not tested!
    if (main.rRec.its.hasSubsurface() &&
        (main.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
      out_veryDirect +=
          main.throughput *
          main.rRec.its.LoSub(scene, main.rRec.sampler, -main.ray.d, 0);
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

    // Some optimizations can be made if this is the last traced segment.
    bool lastSegment = (main.rRec.depth + 1 == m_config->m_maxDepth);

    /* ==================================================================== */
    /*                     Direct illumination sampling                     */
    /* ==================================================================== */

    // Sample incoming radiance from lights (next event estimation).

    {
      const BSDF *mainBSDF = main.rRec.its.getBSDF(main.ray);
      if (main.rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
          mainBSDF->getType() & BSDF::ESmooth &&
          main.rRec.depth + 1 >= m_config->m_minDepth) {

        // Sample an emitter, estimate the direct incident light!

        DirectSamplingRecord dRec(main.rRec.its);
        Point2               lightSample = main.rRec.nextSample2D();
        auto [mainLi, visible] =
            m_scene->sampleEmitterDirectVisible(dRec, lightSample);
        mainLi *= dRec.pdf;

        const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

        BSDFSamplingRecord mainBRec(main.rRec.its,
                                    main.rRec.its.toLocal(dRec.d), ERadiance);
        Spectrum           mainBSDFValue = mainBSDF->eval(mainBRec);
        Float              mainBSDFPdf =
            (emitter->isOnSurface() && dRec.measure == ESolidAngle && visible)
                             ? mainBSDF->pdf(mainBRec)
                             : .0f;

        Float mainDistanceSquared = (main.rRec.its.p - dRec.p).lengthSquared();
        Float mainOpposingCosine  = dot(dRec.n, dRec.d); //? different

        Float mainWeightNumerator = main.pdf * dRec.pdf;
        Float mainWeightDenominator =
            (main.pdf * main.pdf) *
            ((dRec.pdf * dRec.pdf) + (mainBSDFPdf * mainBSDFPdf));

        // ShadowRay check
        if (!m_config->m_strictNormals ||
            dot(main.rRec.its.geoFrame.n, dRec.d) *
                    Frame::cosTheta(mainBRec.wo) >
                .0f) {
          for (int i = 0; i < secondaryCount; ++i) {
            RayState &shifted = shiftedRays[i];

            Spectrum mainContribution(.0f);
            Spectrum shiftedContribution(.0f);
            Float    weight = .0f;

            bool shiftSuccessful = shifted.alive;
            if (shiftSuccessful) {
              // A possible complete offset path

              if (shifted.connection_status ==
                  RAY_CONNECTED /**Already connected */) {
                Float    shiftedBSDFPdf   = mainBSDFPdf;
                Float    shiftedDRecPdf   = dRec.pdf;
                Spectrum shiftedBSDFValue = mainBSDFValue;
                Spectrum shiftedLi        = mainLi;
                Float    jacobian         = 1.f;

                // MIS for the offset path which could also be generated
                // by path tracer in offset pixel

                Float shiftedWeightDenominator =
                    (jacobian * shifted.pdf) * (jacobian * shifted.pdf) *
                    (shiftedDRecPdf * shiftedDRecPdf +
                     shiftedBSDFPdf * shiftedBSDFPdf);
                weight = mainWeightNumerator /
                         (D_EPSILON + shiftedWeightDenominator +
                          mainWeightDenominator);
                mainContribution    = main.throughput * mainBSDFValue * mainLi;
                shiftedContribution = shifted.throughput * shiftedBSDFValue *
                                      shiftedLi * jacobian;
              } else if (shifted.connection_status ==
                         RAY_RECENTLY_CONNECTED /**connected at here */) {
                Vector shiftedWi =
                    normalize(shifted.rRec.its.p - main.rRec.its.p);
                BSDFSamplingRecord shiftBRec(
                    main.rRec.its, main.rRec.its.toLocal(shiftedWi),
                    main.rRec.its.toLocal(dRec.d), ERadiance);
                Float shiftedBSDFPdf =
                    (emitter->isOnSurface() && dRec.measure == ESolidAngle &&
                     visible)
                        ? mainBSDF->pdf(shiftBRec)
                        : 0; // The BSDF sampler can not sample occluded path
                             // segments.

                Float    shiftedDRecPdf   = dRec.pdf;
                Spectrum shiftedBSDFValue = mainBSDF->eval(shiftBRec);
                Spectrum shiftedLi        = mainLi;
                Float    jacobian         = 1.f;

                Float shiftedWeightDenominator =
                    (jacobian * shifted.pdf) * (jacobian * shifted.pdf) *
                    ((shiftedDRecPdf * shiftedDRecPdf) +
                     (shiftedBSDFPdf * shiftedBSDFPdf));
                weight = mainWeightNumerator /
                         (D_EPSILON + shiftedWeightDenominator +
                          mainWeightDenominator);

                mainContribution = main.throughput * (mainBSDFValue * mainLi);
                shiftedContribution = jacobian * shifted.throughput *
                                      (shiftedBSDFValue * shiftedLi);
              } else if (shifted.connection_status ==
                         RAY_NOT_CONNECTED /** no connect */) {

                const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);

                //! The last vertex of the possible offset path is generated by
                //! light sampling, thus only reconnection shift couble be
                //! performed here i.e. only diffuse bsdf would success

                bool       mainAtPointLight = dRec.measure == EDiscrete;
                VertexType mainVertexType =
                               getVertexType(main, *m_config, BSDF::ESmooth),
                           shiftedVertexType =
                               getVertexType(shifted, *m_config, BSDF::ESmooth);
                if (mainAtPointLight ||
                    (mainVertexType == VERTEX_TYPE_DIFFUSE &&
                     shiftedVertexType == VERTEX_TYPE_DIFFUSE)) {

                  DirectSamplingRecord shiftedDRec(shifted.rRec.its);
                  //* same sample2d here, so actually same light sample
                  auto [shiftedLi, shiftedVisible] =
                      m_scene->sampleEmitterDirectVisible(shiftedDRec,
                                                          lightSample);
                  shiftedLi *= shiftedDRec.pdf;
                  Float shiftedDRecPdf = shiftedDRec.pdf;

                  Float shiftedDistanceSquared =
                      (dRec.p - shifted.rRec.its.p).lengthSquared();
                  Vector shiftedWo = normalize(dRec.p - shifted.rRec.its.p);
                  Float  shiftedOpposingCosine = -dot(dRec.n, shiftedWo);

                  BSDFSamplingRecord shiftedBRec{
                      shifted.rRec.its, shifted.rRec.its.toLocal(shiftedWo),
                      ERadiance};
                  // Strict normals check, to make the output match with
                  // bidirectional methods when normal maps are present.
                  if (m_config->m_strictNormals &&
                      dot(shifted.rRec.its.geoFrame.n, shiftedWo) *
                              Frame::cosTheta(shiftedBRec.wo) <
                          0) {
                    // Invalid, non-samplable offset path.
                    shiftSuccessful = false;
                  } else {
                    Spectrum shiftedBSDFValue = shiftedBSDF->eval(shiftedBRec);
                    Float    shiftedBSDFPdf =
                        (emitter->isOnSurface() &&
                         dRec.measure == ESolidAngle && shiftedVisible)
                               ? shiftedBSDF->pdf(shiftedBRec)
                               : 0;
                    Float jacobian =
                        std::abs(shiftedOpposingCosine * mainDistanceSquared) /
                        (D_EPSILON +
                         std::abs(mainOpposingCosine * shiftedDistanceSquared));

                    Float shiftedWeightDenominator =
                        (jacobian * shifted.pdf) * (jacobian * shifted.pdf) *
                        ((shiftedDRecPdf * shiftedDRecPdf) +
                         (shiftedBSDFPdf * shiftedBSDFPdf));
                    weight = mainWeightNumerator /
                             (D_EPSILON + shiftedWeightDenominator +
                              mainWeightDenominator);

                    mainContribution =
                        main.throughput * (mainBSDFValue * mainLi);
                    shiftedContribution = jacobian * shifted.throughput *
                                          (shiftedBSDFValue * shiftedLi);
                  }
                }
              } else {
                std::cout << "Shouldn't arrive here\n";
                exit(1);
              }
            }

            if (!shiftSuccessful) {
              Float shiftedWeightDenominator = .0f;
              weight =
                  mainWeightNumerator / (D_EPSILON + mainWeightDenominator);
              mainContribution    = main.throughput * mainBSDFValue * mainLi;
              shiftedContribution = Spectrum(.0f);
            }

            main.addRadiance(mainContribution, weight);
            shifted.addRadiance(shiftedContribution, weight);
            shifted.addGradient(shiftedContribution - mainContribution, weight);
          } // for loop
        }   // strict normal check
      }
    } // Nee

    /* ==================================================================== */
    /*               BSDF sampling and emitter hits                         */
    /* ==================================================================== */

    BSDFSampleResult mainBSDFResult = sampleBSDF(main);

    if (mainBSDFResult.pdf <= .0f) break;
    const Vector mainWo = main.rRec.its.toWorld(mainBSDFResult.bRec.wo);

    // strict normal check
    if (m_config->m_strictNormals &&
        dot(main.rRec.its.geoFrame.n, mainWo) *
                Frame::cosTheta(mainBSDFResult.bRec.wo) <=
            0) {
      break;
    }

    Intersection previousMainIts = main.rRec.its;
    bool         mainHitEmitter  = false;
    Spectrum     mainLi          = Spectrum(.0f);

    DirectSamplingRecord mainDRec(main.rRec.its);
    const BSDF          *mainBSDF = main.rRec.its.getBSDF(main.ray);

    VertexType mainVertexType =
        getVertexType(main, *m_config, mainBSDFResult.bRec.sampledType);
    VertexType mainNextVertexType;

    main.ray = Ray{main.rRec.its.p, mainWo, main.ray.time};

    // Check if hit the emitter
    if (scene->rayIntersect(main.ray, main.rRec.its)) {
      if (main.rRec.its.isEmitter()) {
        mainLi = main.rRec.its.Le(-main.ray.d);
        mainDRec.setQuery(main.ray, main.rRec.its);
        mainHitEmitter = true;
      }

      // TODO Ignore sub-surface scattering now

      mainNextVertexType =
          getVertexType(main, *m_config, mainBSDFResult.bRec.sampledType);
    } else {
      const Emitter *env = scene->getEnvironmentEmitter();

      if (env) {
        mainLi = env->evalEnvironment(main.ray);
        if (!env->fillDirectSamplingRecord(mainDRec, main.ray)) break;
        mainHitEmitter     = true;
        mainNextVertexType = VERTEX_TYPE_DIFFUSE;
      } else {
        // break
        break;
      }
    }

    // Continue offet path by shift map
    Float mainBSDFPdf     = mainBSDFResult.pdf;
    Float mainPreviousPdf = main.pdf;

    main.throughput *= mainBSDFResult.weight * mainBSDFResult.pdf;
    main.pdf *= mainBSDFResult.pdf;
    main.eta *= mainBSDFResult.bRec.eta;

    const Float mainLumPdf =
        (mainHitEmitter && main.rRec.depth + 1 >= m_config->m_minDepth &&
         !(mainBSDFResult.bRec.sampledType & BSDF::EDelta))
            ? scene->pdfEmitterDirect(mainDRec)
            : .0f;

    Float mainWeightNumerator = mainPreviousPdf * mainBSDFPdf;
    Float mainWeightDenominator =
        (mainPreviousPdf * mainPreviousPdf) *
        ((mainLumPdf * mainLumPdf) + (mainBSDFPdf * mainBSDFPdf));

    for (int i = 0; i < secondaryCount; ++i) {
      RayState &shifted = shiftedRays[i];

      Spectrum shiftedLi{.0f};
      Spectrum mainContribution{.0f};
      Spectrum shifedContribution{.0f};
      Float    weight = .0f;

      bool postponedShiftEnd = false;

      if (shifted.alive) {
        Float shiftedPreviousPdf = shifted.pdf;

        if (shifted.connection_status == RAY_CONNECTED /** connected */) {
          // TODO
        } else if (shifted.connection_status ==
                   RAY_RECENTLY_CONNECTED /** connect here */) {
          // TODO
        } else if (shifted.connection_status ==
                   RAY_NOT_CONNECTED /** shift map */) {
          // TODO

          // Random number replay shitmap
          // TODO
        } else {
          std::cout << "Shouldn't arrive here\n";
          exit(1);
        }
      }

      // TODO handle failed shifted path
      // TODO account the contribution (possible)
    }

    // Terminate the path if not intersect the scene
    // Apply rr
  }
}

MTS_NAMESPACE_END