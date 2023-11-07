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

    Spectrum L(.0f), beta(1.f), r_u(1.f), r_l(1.f);
    Float    prev_scatter_pdf = .0f;
    Point    prev_p;
    Normal   prev_n;

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

        Spectrum T_maj = SampleTrMajorant(medium, ray, tMax, rng,
                                          [&](MajorantSamplingRecord majRec) {
                                            //
                                            return false;
                                          });
      }
    }
  }
};

MTS_NAMESPACE_END
