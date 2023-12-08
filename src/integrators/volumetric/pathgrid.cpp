// clang-format off

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

#include "pathgrid/pathgrid.h"
// clang-format on

/**
 * A volumetric path-tracer using subpath reuse for multiple scattering
 */

MTS_NAMESPACE_BEGIN

constexpr uint32_t K = 4;

class PathGrid : public MonteCarloIntegrator {
public:
  PathGrid(const Properties &props) : MonteCarloIntegrator(props) {
    //
  }

  PathGrid(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {
    Log(EError, "PathGrid serialization is not support\n");
  }

  void serialize(Stream *stream, InstanceManager *manager) const override {
    MonteCarloIntegrator::serialize(stream, manager);
    Log(EError, "Path serialization is not support\n");
  }

  virtual Spectrum Li(const RayDifferential &r,
                      RadianceQueryRecord   &rRec) const override {
    using RNG = RNG_CX;

    // Uniformly sample a hero wave length
    const uint32_t heroChannel = SampleHeroChannel(rRec.nextSample1D());

    const Scene    *scene = rRec.scene;
    Intersection   &its   = rRec.its;
    RayDifferential ray(r);

    Spectrum L(.0f), beta(1.f), r_u(1.f), r_l(1.f);

    Float  prev_scatter_pdf = .0f;
    Point  prev_p;
    Normal prev_n;
    bool   specular_bounce = true;

    while (true) {
      if (beta.isZero()) break;

      ray.heroChannel = heroChannel;
      scene->rayIntersect(ray, its);

      if (const Medium *medium = rRec.medium; medium) {
        bool scattered  = false;
        bool terminated = false;

        Float    t_max = its.isValid() ? its.t : INFINITY;
        uint64_t seed  = ((uint64_t)rand() << 32) + (uint64_t)rand();
        RNG      rng(seed);

        Spectrum T_maj = SampleTrMajorant(
            medium, ray, t_max, rng, ETrackingType::ENaiveDeltaTracking,
            [&](MajorantSamplingRecord majRec) {
              // Sample a type of collision
              Spectrum sigma_a = majRec.sigma_a;
              Spectrum sigma_s = majRec.sigma_s;
              Spectrum sigma_n = majRec.sigma_n;
              Spectrum tr_maj  = majRec.tr_majorant;

              ECollisionType scatter_type = SampleCollision(
                  sigma_a, sigma_s, sigma_n, rng.randomFloat(), heroChannel);

              if (scatter_type & ECollisionType::Absorb) {
                // TODO Add volumetric emission
                terminated = true;
                return false;
              } else if (scatter_type & ECollisionType::Scatter) {
                if (rRec.depth++ >= m_maxDepth) {
                  terminated = true;

                  // TODO multiple-scattering should be added here
                  uint32_t idxs[K];

                  return false;
                }

                // update the throughput
                Float pdf = tr_maj[heroChannel] * sigma_s[heroChannel];
                beta *= tr_maj * sigma_s / pdf;
                r_u *= tr_maj * sigma_s / pdf;

                DirectSamplingRecord dRec{majRec.p, .0f};
                dRec.refN = Normal(.0f);

                L += beta * SampleVolumetricNEE(scene, dRec, -ray.d, medium,
                                                rRec.sampler, r_u, heroChannel,
                                                medium->getPhaseFunction(),
                                                nullptr);
                // sample a new direction
                const PhaseFunction *phase = medium->getPhaseFunction();
                MediumSamplingRecord
                    mRec; // TODO mRec will be used in some phase
                PhaseFunctionSamplingRecord pRec{mRec, -ray.d};
                Float                       phasePdf;
                Float phaseWeight = phase->sample(pRec, phasePdf, rRec.sampler);

                if (phaseWeight == .0f || phasePdf == .0f)
                  terminated = true;
                else {
                  scattered       = true;
                  specular_bounce = false;

                  beta     = beta * phaseWeight;
                  r_l      = r_u / phasePdf;
                  ray      = Ray{majRec.p, pRec.wo, .0f};
                  ray.mint = Epsilon;

                  prev_p = majRec.p;
                  prev_n = Normal(.0f);
                }
                return false;
              } else if (scatter_type & ECollisionType::Null) {
                //* A null scatter occurs, just keep tracking
                Spectrum sigma_maj = majRec.sigma_maj;

                Float pdf = tr_maj[heroChannel] * sigma_n[heroChannel];
                beta *= tr_maj * sigma_n / pdf;

                r_u *= tr_maj * sigma_n / pdf;
                r_l *= tr_maj * sigma_maj / pdf;

                ray      = Ray{majRec.p, ray.d, .0f};
                ray.mint = Epsilon;

                return !beta.isZero() && !r_u.isZero();
              } else {
                Log(EError, "Shouldn't arrive here");
                exit(1);
              }
            });

        if (terminated) break;
        if (scattered) continue;

        beta *= T_maj / T_maj[heroChannel];
        r_u *= T_maj / T_maj[heroChannel];
        r_l *= T_maj / T_maj[heroChannel];
      }

      if (!its.isValid()) {
        Spectrum Le = scene->evalEnvironment(ray);
        if (rRec.depth == 1 || specular_bounce) {
          // First hit or pure specular bounce
          L += beta * Le / r_u.average();
        } else {
          // Apply MIS
          Float p_l = PDF_Nee(scene, prev_p, prev_n, ray);
          r_l       = r_l * p_l;
          L += beta * Le / (r_u + r_l).average();
        }
        break;
      }

      if (const Emitter *emitter = its.shape->getEmitter(); emitter) {
        Spectrum Le = emitter->eval(its, -ray.d);
        if (rRec.depth == 1 || specular_bounce) {
          // First hit or pure specular bounce
          L += beta * Le / r_u.average();
        } else {
          // Apply MIS
          Float p_l = PDF_Nee(scene, prev_p, prev_n, ray, &its);
          r_l       = r_l * p_l;
          L += beta * Le / (r_u + r_l).average();
        }
      }

      const BSDF *bsdf = its.getBSDF(ray);
      if (bsdf->isNull()) {
        ray      = Ray{its.p + Epsilon * ray.d, ray.d, ray.time};
        ray.mint = Epsilon;
        if (its.isMediumTransition()) rRec.medium = its.getTargetMedium(ray.d);
        continue;
      }

      if (rRec.depth++ >= m_maxDepth) break;

      DirectSamplingRecord dRec(its);
      if (bsdf->getType() & BSDF::ESmooth)
        L += beta * SampleVolumetricNEE(scene, dRec, -ray.d, rRec.medium,
                                        rRec.sampler, r_u, heroChannel, nullptr,
                                        bsdf, &its);

      // BSDF sampling
      BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
      Float              bsdfPdf;
      Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
      if (bsdfWeight.isZero()) break;

      const Vector wo        = its.toWorld(bRec.wo);
      Float        woDotGeoN = dot(its.geoFrame.n, wo);
      if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals) break;

      beta *= bsdfWeight;
      ray = Ray(its.p + Epsilon * wo, wo, Epsilon, INFINITY, ray.time);
      specular_bounce  = (bRec.sampledType &
                         (BSDF::EDeltaReflection | BSDF::EDeltaTransmission));
      r_l              = r_u / bsdfPdf;
      prev_scatter_pdf = bsdfPdf;

      prev_p = its.p;
      prev_n = its.shFrame.n;

      if (its.isMediumTransition()) rRec.medium = its.getTargetMedium(ray.d);

      if (rRec.depth >= m_rrDepth) {
        Float q = std::min(beta.max(), (Float)0.95f);
        if (rRec.nextSample1D() >= q) break;
        beta /= q;
      }
    }

    return L;
  }

protected:
  enum ECollisionType { Absorb = 0b0001, Scatter = 0b0010, Null = 0b0100 };

  ECollisionType SampleCollision(const Spectrum &sigma_a,
                                 const Spectrum &sigma_s,
                                 const Spectrum &sigma_n, Float u,
                                 uint32_t channel) const {
    // Only first channel here
    Float sum = sigma_a[channel] + sigma_s[channel] + sigma_n[channel];
    Float P_a = sigma_a[channel] / sum;
    Float P_s = sigma_s[channel] / sum;

    if (u < P_a) return ECollisionType::Absorb;
    if (u < P_a + P_s) return ECollisionType::Scatter;
    return ECollisionType::Null;
  }

  uint32_t SampleHeroChannel(Float u) const {
    return math::clamp((int)(u * 3), 0, 2);
  }

  Spectrum SampleVolumetricNEE(const Scene *scene, DirectSamplingRecord &dRec,
                               const Vector &wi, const Medium *medium,
                               Sampler *sampler, Spectrum r_p,
                               const uint32_t       heroChannel,
                               const PhaseFunction *phase = nullptr,
                               const BSDF          *bsdf  = nullptr,
                               const Intersection  *its   = nullptr) const {
    using RNG = RNG_CX;

    Spectrum Le{.0f};
    Spectrum scatterVal{.0f};
    Float    scatterPdf = .0f;
    Point2   sample     = sampler->next2D();

    Le = scene->sampleEmitterDirect(dRec, sample, false);

    if (Le.isZero() || dRec.pdf == .0f) return Spectrum(.0f);

    const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

    if (bsdf && (bsdf->getType() & BSDF::ESmooth)) {
      BSDFSamplingRecord bRec(*its, its->toLocal(dRec.d), ERadiance);
      scatterVal = bsdf->eval(bRec);
      scatterPdf = bsdf->pdf(bRec);
    } else if (phase) {
      MediumSamplingRecord        mRec; // TODO
      PhaseFunctionSamplingRecord pRec{mRec, wi, dRec.d};
      scatterVal = Spectrum{phase->eval(pRec)};
      scatterPdf = phase->pdf(pRec);
    }

    Spectrum        tr(1.f), r_l(1.f), r_u(1.f);
    RayDifferential shadowRay{dRec.ref, dRec.d, .0f};
    shadowRay.mint        = Epsilon;
    shadowRay.maxt        = dRec.dist * (1 - ShadowEpsilon);
    shadowRay.heroChannel = heroChannel;

    uint64_t      seed = ((uint64_t)rand() << 32) + (uint64_t)rand();
    RNG           rng{seed};
    const Medium *currentMedium = medium;
    Intersection  si;
    while (true) {

      scene->rayIntersect(shadowRay, si);

      if (si.isValid()) {
        if (const BSDF *bsdf = si.getBSDF(); !bsdf->isNull()) {
          // Hit opaque surface
          tr = Spectrum(.0f);
          break;
        }
      }

      if (currentMedium) {
        // Accumulate tr
        Float    tMax = si.isValid() ? si.t : shadowRay.maxt;
        Spectrum T_maj =
            SampleTrMajorant(currentMedium, shadowRay, tMax, rng,
                             ETrackingType::ENaiveDeltaTracking,
                             [&](MajorantSamplingRecord majRec) {
                               Spectrum sigma_n   = majRec.sigma_n;
                               Spectrum sigma_maj = majRec.sigma_maj;
                               Spectrum tr_maj    = majRec.tr_majorant;
                               Float    pdf =
                                   tr_maj[heroChannel] * sigma_maj[heroChannel];

                               tr *= tr_maj * sigma_n / pdf;
                               r_l *= tr_maj * sigma_maj / pdf;
                               r_u *= tr_maj * sigma_n / pdf;

                               if (tr.isZero()) return false;
                               return true;
                             });

        tr *= T_maj / T_maj[heroChannel];
        r_l *= T_maj / T_maj[heroChannel];
        r_u *= T_maj / T_maj[heroChannel];
      }

      if (tr.isZero()) break;

      shadowRay.o    = si.p;
      shadowRay.mint = Epsilon;
      shadowRay.maxt = shadowRay.maxt - si.t;

      if (!si.isValid()) break;

      if (si.isMediumTransition())
        currentMedium = si.getTargetMedium(shadowRay.d);
    }

    r_l *= r_p;
    r_u *= r_p * scatterPdf / dRec.pdf;

    if (!(emitter->isOnSurface() && dRec.measure == ESolidAngle))
      return scatterVal * tr * Le / r_l.average();
    else
      return scatterVal * tr * Le / (r_l + r_u).average();
  }

  Float PDF_Nee(const Scene *scene, Point prev_p, Normal prev_n,
                const RayDifferential &ray,
                const Intersection    *its = nullptr) const {
    DirectSamplingRecord dRec;
    Float                pdf = .0f;

    if (its) {
      dRec.ref  = prev_p;
      dRec.refN = prev_n;

      if (dot(dRec.refN, ray.d) < 0) dRec.refN *= -1;

      dRec.p       = its->p;
      dRec.n       = its->shFrame.n;
      dRec.measure = ESolidAngle;
      dRec.uv      = its->uv;
      dRec.object  = its->shape->getEmitter();
      dRec.d       = normalize(dRec.p - prev_p);
      dRec.dist    = (dRec.p - prev_p).length();

      pdf = scene->pdfEmitterDirect(dRec);
    } else {
      const Emitter *env = scene->getEnvironmentEmitter();
      if (env && env->fillDirectSamplingRecord(dRec, ray))
        return scene->pdfEmitterDirect(dRec);
    }

    return pdf;
  }

private:
  std::unique_ptr<pathgrid::PathGrid> m_pathGrid;

  MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS(PathGrid, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(PathGrid, "Volumetric PathGrid")

MTS_NAMESPACE_END