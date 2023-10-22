#include <mitsuba/core/statistics.h>
#include <mitsuba/render/scene.h>

/**
    Volumetric path tracer based on
    "A null-scattering path integral formulation of light transport"
*/

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Volumetric null-scatter path tracer", "Average path length",
                                  EAverage);

struct ScatterRecord {};

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

    int bounces = 0;

    while (true) {
      bool found_intersection = rRec.rayIntersect(ray);

      if (rRec.medium) {
        bool  scattered  = false;
        bool  terminated = false;
        Float tmax       = found_intersection ? its.t : INFINITY;

        RNG_CX rng;

        SampleTrMajorant(rRec.medium, ray, tmax, rng, [&](MajorantSamplingRecord maj_rec) {
          Spectrum tr_maj    = maj_rec.tr_majorant;
          Spectrum sigma_maj = maj_rec.sigma_maj;

          beta *= tr_maj * sigma_maj / (tr_maj[0] * sigma_maj[0]);

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

            //* Luminaire sampling

            if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
              L += beta * SampleLd(maj_rec.p, scene, maj_rec.medium, rRec.nextSample2D(), rng);
            }

          }

          else /*Null-collision*/ {
            return true;
          }
        });
      }
    }

    return L;
  }
  MTS_DECLARE_CLASS()

protected:
  Spectrum SampleLd(Point3 p, const Scene *scene, const Medium *medium, Point2 uLight,
                    const RNG_CX &rng) const {
    DirectSamplingRecord dRec(p, .0f);
    Spectrum             Ld = scene->SampleVolumetricDirect(&dRec, medium, uLight, rng);

    if (!Ld.isZero()) {
      const Emitter *emitter = static_cast<const Emitter *>(dRec.object);
    }
  }
};

MTS_IMPLEMENT_CLASS_S(VolumetricNSPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(VolumetricNSPathTracer, "Volumetric null-scattering path tracer");

MTS_NAMESPACE_END