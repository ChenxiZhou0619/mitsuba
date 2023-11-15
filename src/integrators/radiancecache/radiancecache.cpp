// clang-format off
#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

// clang-format on

MTS_NAMESPACE_BEGIN

// KD-tree + vMF radiance cache

class RadianceCachePathTracer : public MonteCarloIntegrator {
public:
  RadianceCachePathTracer(const Properties &props)
      : MonteCarloIntegrator(props) {
    //
  }
};

MTS_NAMESPACE_END