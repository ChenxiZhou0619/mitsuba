// clang-format off

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

// clang-format on

MTS_NAMESPACE_BEGIN

class SDTrackingPT : public MonteCarloIntegrator {
public:
  SDTrackingPT(const Properties &props) : MonteCarloIntegrator(props) {
    //
  }

  SDTrackingPT(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {
    Log(EError, "SDTrackingPT serialization is not support\n");
  }

  void serialize(Stream *stream, InstanceManager *manager) const override {
    MonteCarloIntegrator::serialize(stream, manager);
    Log(EError, "SDTrackingPT serialization is not support\n");
  }

  virtual Spectrum Li(const RayDifferential &r,
                      RadianceQueryRecord   &rRec) const override {
    using RNG = RNG_CX;

    const Scene    *scene = rRec.scene;
    Intersection   &its   = rRec.its;
    RayDifferential ray(r);

    Spectrum        L(.0f), beta(1.f);
    constexpr Float r_u = 1.f;
    Float           r_l = 1.f;

    Point  prev_p;
    Normal prev_n;

    bool specular_bounce = true;

    while (true) {
      if (beta.isZero()) break;

      scene->rayIntersect(ray, its);

      if (const Medium *medium = rRec.medium; medium) {
        bool scattered  = false;
        bool terminated = false;

        Float tMax = its.isValid() ? its.t : INFINITY;

        uint64_t seed = ((uint64_t)rand() << 32) + (uint64_t)rand();
        RNG      rng(seed);

        Float    pdfPassThrough;
        Spectrum T_maj = SampleTrMajorant(
            medium, ray, tMax, rng, ETrackingType::ESpectralTracking,
            [&](MajorantSamplingRecord majRec) {
              Spectrum sigma_a   = majRec.sigma_a;
              Spectrum sigma_s   = majRec.sigma_s;
              Spectrum sigma_n   = majRec.sigma_n;
              Spectrum tr_maj    = majRec.tr_majorant;
              Float    pdfFlight = majRec.pdf_flight;

              Float          pdfCollision;
              ECollisionType scatterType =
                  SampleCollision(sigma_a, sigma_s, sigma_n, beta,
                                  rng.randomFloat(), &pdfCollision);

              if (scatterType & ECollisionType::Absorb) {
                terminated = true;
                return false;
              } else if (scatterType & ECollisionType::Scatter) {
                // Sample a real scatter

                if (rRec.depth++ >= m_maxDepth) {
                  terminated = true;
                  return false;
                }

                // update the throughput
                beta *= tr_maj * sigma_s / (pdfFlight * pdfCollision);

                if (!beta.isValid()) {
                  printf("Stop here 1\n");
                }

                DirectSamplingRecord dRec{majRec.p, .0f};
                dRec.refN = Normal(.0f);

                L += beta * SampleVolumetricNEE(scene, dRec, -ray.d, medium,
                                                rRec.sampler, r_u,
                                                medium->getPhaseFunction());

                const PhaseFunction        *phase = medium->getPhaseFunction();
                MediumSamplingRecord        mRec; // TODO mRec
                PhaseFunctionSamplingRecord pRec{mRec, -ray.d};

                Float phasePdf;
                Float phaseWeight = phase->sample(pRec, phasePdf, rRec.sampler);

                if (phaseWeight == .0f || phasePdf == .0f) {
                  terminated = true;
                } else {
                  scattered       = true;
                  specular_bounce = false;

                  beta     = beta * phaseWeight;
                  r_l      = r_u / phasePdf;
                  ray.o    = majRec.p;
                  ray.d    = pRec.wo;
                  ray.mint = Epsilon;

                  prev_p = majRec.p;
                  prev_n = Normal(.0f);
                }
                return false;
              } else if (scatterType & ECollisionType::Null) {
                // Sample a null collision

                // update the throughput
                beta *= tr_maj * sigma_n / (pdfFlight * pdfCollision);

                r_l *= (1.f / pdfCollision);

                if (!beta.isValid()) {
                  printf("Stop here 2\n");
                }

                ray      = Ray{majRec.p, ray.d, .0f};
                ray.mint = Epsilon;

                return !beta.isZero();
              } else {
                Log(EError, "Shouldn't arrive here");
                exit(1);
              }
            },
            &pdfPassThrough);

        if (terminated) break;
        if (scattered) continue;

        beta *= T_maj / pdfPassThrough;
      }

      if (!its.isValid()) {
        Spectrum Le = scene->evalEnvironment(ray);
        if (rRec.depth == 1 || specular_bounce) {
          // First hit or pure specular bounce
          L += beta * Le / r_u;
        } else {
          // Apply MIS
          Float p_l = PDF_Nee(scene, prev_p, prev_n, ray);
          r_l       = r_l * p_l;
          L += beta * Le / (r_u + r_l);
        }
        break;
      }

      if (const Emitter *emitter = its.shape->getEmitter(); emitter) {
        Spectrum Le = emitter->eval(its, -ray.d);
        if (rRec.depth == 1 || specular_bounce) {
          // First hit or pure specular bounce
          L += beta * Le / r_u;
        } else {
          // Apply MIS
          Float p_l = PDF_Nee(scene, prev_p, prev_n, ray, &its);
          r_l       = r_l * p_l;
          L += beta * Le / (r_u + r_l);
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
                                        rRec.sampler, r_u, nullptr, bsdf, &its);

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
      specular_bounce = (bRec.sampledType &
                         (BSDF::EDeltaReflection | BSDF::EDeltaTransmission));
      r_l             = r_u / bsdfPdf;

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
                                 const Spectrum &sigma_n, const Spectrum &beta,
                                 Float u, Float *pdf) const {
     //Float sigAMax = sigma_a.max();
     //Float sigSMax = sigma_s.max();
     //Float sigNMax = sigma_n.max();
    //
     //Float invC = 1.f / (sigAMax + sigSMax + sigNMax);
    //
     //Float Pa = sigAMax * invC;
     //Float Ps = sigSMax * invC;
     //Float Pn = sigNMax * invC;

    Float Ca = (sigma_a * beta).average();
    Float Cs = (sigma_s * beta).average();
    Float Cn = (sigma_n * beta).average();

    Float invC = 1.f / (Ca + Cs + Cn);
    Float Pa   = Ca * invC;
    Float Ps   = Cs * invC;
    Float Pn   = Cn * invC;

    if (u < Pa) {
      *pdf = Pa;
      return ECollisionType::Absorb;
    } else if (u < Pa + Ps) {
      *pdf = Ps;
      return ECollisionType::Scatter;
    } else {
      *pdf = Pn;
      return ECollisionType::Null;
    }
  }

  Spectrum SampleVolumetricNEE(const Scene *scene, DirectSamplingRecord &dRec,
                               const Vector &wi, const Medium *medium,
                               Sampler *sampler, Float r_p,
                               const PhaseFunction *phase = nullptr,
                               const BSDF          *bsdf  = nullptr,
                               const Intersection  *its   = nullptr) const {
    using RNG = RNG_CX;

    Spectrum Le{.0f};
    Spectrum scatterVal{.0f};
    Float    scatterPdf = .0f;
    Point2   sample     = sampler->next2D();

    Le = scene->sampleEmitterDirect(dRec, sample, false);

    if (Le.isZero() || dRec.pdf == .0f) return Spectrum{.0f};

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

    Spectrum tr(1.f);
    Float    r_u = 1.f;
    Float    r_l = 1.f;

    RayDifferential shadowRay{dRec.ref, dRec.d, .0f};
    shadowRay.mint = Epsilon;
    shadowRay.maxt = dRec.dist * (1 - ShadowEpsilon);

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
        Float    pdfPassThrough;
        Spectrum T_maj = SampleTrMajorant(
            currentMedium, shadowRay, tMax, rng,
            ETrackingType::ESpectralTracking,
            [&](MajorantSamplingRecord majRec) {
              Float sigAMax = majRec.sigma_a.max();
              Float sigSMax = majRec.sigma_s.max();
              Float sigNMax = majRec.sigma_n.max();

              Float invC = 1.f / (sigAMax + sigSMax + sigNMax);

              Float Pn = sigNMax * invC;

              tr *= majRec.tr_majorant * majRec.sigma_n / majRec.pdf_flight;
              r_u *= Pn;

              if (!tr.isValid()) {
                printf("Stop here 4\n");
              }

              if (tr.isZero()) return false;
              return true;
            },
            &pdfPassThrough);

        if (tr.isZero()) break;

        tr *= T_maj / pdfPassThrough;
        if (!tr.isValid()) {
          printf("Stop here 5\n");
        }
      }

      shadowRay.o    = si.p;
      shadowRay.mint = Epsilon;
      shadowRay.maxt = shadowRay.maxt - si.t;

      if (!si.isValid()) break;

      if (si.isMediumTransition())
        currentMedium = si.getTargetMedium(shadowRay.d);
    }

    r_l *= r_p;
    r_u *= r_p * scatterPdf / dRec.pdf;

    if (!tr.isValid()) {
      printf("Stop here 3\n");
    }

    if (!(emitter->isOnSurface() && dRec.measure == ESolidAngle))
      return scatterVal * tr * Le / r_l;
    else
      return scatterVal * tr * Le / (r_l + r_u);
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
        pdf = scene->pdfEmitterDirect(dRec);
    }
    return pdf;
  }

  MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS(SDTrackingPT, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(SDTrackingPT, "Spectral tracking path tracer")

MTS_NAMESPACE_END
