#include <mitsuba/core/statistics.h>
#include <mitsuba/render/scene.h>
#include <variant>
/**
    Volumetric path tracer based on
    "A Null-scattering Path Integral Formulation of Light Transport"
*/

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Volumetric null-scatter path tracer", "Average path length",
                                  EAverage);

struct ScatterRecord {
  std::variant<const BSDF *, const PhaseFunction *> scatterFunction;
  Vector3                                           localWi, localWo;
  Point3                                            p;

  ScatterRecord(const Intersection &its, const BSDF *bsdf);

  ScatterRecord(const MajorantSamplingRecord &mRec, Vector3 wi);
};

class VolumetricNSPathTracer : public MonteCarloIntegrator {
public:
  VolumetricNSPathTracer(const Properties &props) : MonteCarloIntegrator(props){};

  VolumetricNSPathTracer(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {}

  Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
    const Scene  *scene = rRec.scene;
    Intersection &its   = rRec.its;

    RayDifferential ray(r);
    Spectrum        L(.0f);
    Spectrum        beta(1.f);
    Spectrum        r_u(1.f); // P_uni(lambda_0, ..., lambda_n) / P_uni(lambda)
    Spectrum        r_l(1.f);

    int  bounces         = 0;
    bool specular_bounce = false;

    while (true) {
      bool found_intersection = rRec.rayIntersect(ray);

      if (const Medium *medium = rRec.medium; medium) {
        bool  scattered  = false;
        bool  terminated = false;
        Float tmax       = found_intersection ? its.t : INFINITY;

        RNG_CX rng;

        Spectrum tMaj =
            SampleTrMajorant(medium, ray, tmax, rng, [&](MajorantSamplingRecord maj_rec) {
              Spectrum tr_maj    = maj_rec.tr_majorant;
              Spectrum sigma_maj = maj_rec.sigma_maj;
              Spectrum sigma_s   = maj_rec.sigma_s;
              Spectrum sigma_n   = maj_rec.sigma_n;

              //* Sample the collision type
              Float p_absorb  = maj_rec.sigma_a[0] / sigma_maj[0];
              Float p_scatter = maj_rec.sigma_s[0] / sigma_maj[0];

              auto SampleCollision = [&](Float u, Float p_absorb, Float p_scatter) -> int {
                if (u < p_absorb)
                  return 0;
                else if (u < p_absorb + p_scatter)
                  return 1;
                return 2;
              };

              int collision_type = SampleCollision(rng.randomFloat(), p_absorb, p_scatter);

              if (collision_type == 0 /* Absorbtion*/) {
                terminated = true;
                return false; // Return false to terminate the tracking
              }

              else if (collision_type == 1 /* Scattering*/) {
                if (bounces >= m_maxDepth) {
                  terminated = true;
                  return false;
                }

                Float pdf = tr_maj[0] * sigma_s[0];
                beta *= tr_maj * sigma_s / pdf;
                r_u *= tr_maj * sigma_s / pdf;

                //* Luminaire sampling
                ScatterRecord sRec{maj_rec, -ray.d};
                L += beta * SampleLd(sRec, scene, medium, rRec.nextSample2D(), rng, r_u);

                //* Phase sampling
                const PhaseFunction *phase = medium->getPhaseFunction();
                MediumSamplingRecord mRec; // TODO some phase sample may need mRec information
                PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                Float                       phasePdf;
                Float phaseWeight = phase->sample(pRec, phasePdf, rRec.sampler);

                beta *= phaseWeight;
                r_l /= phasePdf;

                scattered = true;
                ray       = Ray(maj_rec.p, pRec.wo, ray.time);
                ray.mint  = .0f;
                return false;
              }

              else /*Null-collision*/ {

                Float pdf = tr_maj[0] * sigma_n[0];
                beta *= tr_maj * sigma_n / pdf;
                r_u *= tr_maj * sigma_n / pdf;
                r_l *= tr_maj * sigma_maj / pdf;

                return true;
              }
            });

        if (terminated || beta.isZero()) break;
        if (scattered) continue;

        //* sample pass through the media
        beta *= tMaj / tMaj[0];
      }

      if (!found_intersection) {
        // TODO add environment light
        break;
      }

      if (const Emitter *emitter = its.shape->getEmitter(); emitter) {
        Spectrum Le = its.Le(-ray.d);

        if (bounces == 0 || specular_bounce) {
          L += beta * Le * r_u.average();
        } else {
        }
      }
    }

    return L;
  }
  MTS_DECLARE_CLASS()

protected:
  Spectrum SampleLd(ScatterRecord s, const Scene *scene, const Medium *medium, Point2 uLight,
                    const RNG_CX &rng, Spectrum r_u, const Intersection *its = nullptr) const {
    DirectSamplingRecord dRec(s.p, .0f);
    Spectrum             Ld = scene->SampleVolumetricDirect(&dRec, medium, uLight, rng);

    if (!Ld.isZero()) {
      const Emitter *emitter = static_cast<const Emitter *>(dRec.object);
      Spectrum       fs(.0f);
      Float          pdfScatter;

      if (s.scatterFunction.index() == 0) {
        const BSDF        *bsdf = std::get<const BSDF *>(s.scatterFunction);
        BSDFSamplingRecord bRec(*its, its->toLocal(dRec.d));
        fs = bsdf->eval(bRec);
        pdfScatter =
            (emitter->isOnSurface() && dRec.measure == ESolidAngle) ? bsdf->pdf(bRec) : .0f;
      } else {
        const PhaseFunction        *phase = std::get<const PhaseFunction *>(s.scatterFunction);
        MediumSamplingRecord        mRec; // TODO some phase eval may need mRec information
        PhaseFunctionSamplingRecord pRec{mRec, s.localWi, dRec.d};
        fs = Spectrum(phase->eval(pRec));
        pdfScatter =
            (emitter->isOnSurface() && dRec.measure == ESolidAngle) ? phase->pdf(pRec) : .0f;
      }

      if (!fs.isZero()) {
        Float misw = dRec.pdf / (dRec.pdf + pdfScatter);
        return misw * fs * Ld / dRec.pdf;
      }
    }
  }
};

MTS_IMPLEMENT_CLASS_S(VolumetricNSPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(VolumetricNSPathTracer, "Volumetric null-scattering path tracer");

MTS_NAMESPACE_END