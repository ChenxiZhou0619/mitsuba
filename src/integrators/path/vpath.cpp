#include <mitsuba/core/statistics.h>
#include <mitsuba/render/scene.h>

/**
    Volumetric path tracer based on
    "A null-scattering path integral formulation of light transport"
*/

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Volumetric null-scatter path tracer", "Average path length",
                                  EAverage);

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

    while (true) {
      bool found_intersection = rRec.rayIntersect(ray);

      if (rRec.medium) {
        //
      }

      if (found_intersection) {
        L.fromLinearRGB(its.geoFrame.n.x, its.geoFrame.n.y, its.geoFrame.n.z);
        L += Spectrum(1.f);
        L *= .5f;
      }
      break;
    }

    return L;
  }

  MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(VolumetricNSPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(VolumetricNSPathTracer, "Volumetric null-scattering path tracer");

MTS_NAMESPACE_END