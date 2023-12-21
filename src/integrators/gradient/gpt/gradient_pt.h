#pragma once

// clang-format off

#include <mitsuba/mitsuba.h>
#include "gpt_wr.h"
#include "gpt_proc.h"

// clang-format on
MTS_NAMESPACE_BEGIN

/**
 * ------------------------- Declerations -------------------------
 */

/// A threshold to use in positive denominators to avoid division by zero.
const Float D_EPSILON = (Float)(1e-14);

/// Classification of vertices into diffuse and glossy.
enum VertexType {
  VERTEX_TYPE_GLOSSY, ///< "Specular" vertex that requires the half-vector
                      ///< duplication shift.
  VERTEX_TYPE_DIFFUSE ///< "Non-specular" vertex that is rough enough for the
                      ///< reconnection shift.
};

enum RayConnection {
  RAY_NOT_CONNECTED,      ///< Not yet connected - shifting in progress.
  RAY_RECENTLY_CONNECTED, ///< Connected, but different incoming direction so
                          ///< needs a BSDF evaluation.
  RAY_CONNECTED ///< Connected, allows using BSDF values from the base path.
};

/// Stores the results of a BSDF sample.
/// Do not confuse with Mitsuba's BSDFSamplingRecord.
struct BSDFSampleResult {
  BSDFSamplingRecord bRec;   ///< The corresponding BSDF sampling record.
  Spectrum           weight; ///< BSDF weight of the sampled direction.
  Float              pdf;    ///< PDF of the BSDF sample.
};

/// Describes the state of a ray that is being traced in the scene.
struct RayState {
  RayState()
      : radiance(0.0f), gradient(0.0f), eta(1.0f), pdf(1.0f),
        throughput(Spectrum(0.0f)), alive(true),
        connection_status(RAY_NOT_CONNECTED) {}

  /// Adds radiance to the ray.
  inline void addRadiance(const Spectrum &contribution, Float weight) {
    Spectrum color = contribution * weight;
    radiance += color;
  }

  /// Adds gradient to the ray.
  inline void addGradient(const Spectrum &contribution, Float weight) {
    Spectrum color = contribution * weight;
    gradient += color;
  }

  RayDifferential ray; ///< Current ray.

  Spectrum throughput; ///< Current throughput of the path.
  Float    pdf;        ///< Current PDF of the path.

  // Note: Instead of storing throughput and pdf, it is possible to store
  // Veach-style weight (throughput divided by pdf), if relative PDF (offset_pdf
  // divided by base_pdf) is also stored. This might be more stable numerically.

  Spectrum radiance; ///< Radiance accumulated so far.
  Spectrum gradient; ///< Gradient accumulated so far.

  RadianceQueryRecord rRec; ///< The radiance query record for this ray.
  Float               eta;  ///< Current refractive index of the ray.
  bool alive; ///< Whether the path matching to the ray is still good. Otherwise
              ///< it's an invalid offset path with zero PDF and throughput.

  RayConnection connection_status; ///< Whether the ray has been connected to
                                   ///< the base path, or is in progress.
};

/**
 * ------------------------- Helper Functions -------------------------
 */

/// Returns the vertex type of a vertex by its roughness value.
inline VertexType
getVertexTypeByRoughness(Float                           roughness,
                         const GradientPathTracerConfig &config) {
  if (roughness <= config.m_shiftThreshold) {
    return VERTEX_TYPE_GLOSSY;
  } else {
    return VERTEX_TYPE_DIFFUSE;
  }
}

/// Returns whether the given ray sees the environment.
inline bool testEnvironmentVisibility(const Scene *scene, const Ray &ray) {
  const Emitter *env = scene->getEnvironmentEmitter();
  if (!env) {
    return false;
  }

  Ray shadowRay(ray);
  shadowRay.setTime(ray.time);
  shadowRay.setOrigin(ray.o);
  shadowRay.setDirection(ray.d);

  DirectSamplingRecord directSamplingRecord;
  env->fillDirectSamplingRecord(directSamplingRecord, shadowRay);

  shadowRay.mint = Epsilon;
  shadowRay.maxt = ((Float)1.0 - ShadowEpsilon) * directSamplingRecord.dist;

  return !scene->rayIntersect(shadowRay);
}

/// Returns whether point1 sees point2.
inline bool testVisibility(const Scene *scene, const Point3 &point1,
                           const Point3 &point2, Float time) {
  Ray shadowRay;
  shadowRay.setTime(time);
  shadowRay.setOrigin(point1);
  shadowRay.setDirection(point2 - point1);
  shadowRay.mint = Epsilon;
  shadowRay.maxt = (Float)1.0 - ShadowEpsilon;

  return !scene->rayIntersect(shadowRay);
}

/// Returns the vertex type (diffuse / glossy) of a vertex, for the purposes of
/// determining the shifting strategy.
///
/// A bare classification by roughness alone is not good for multi-component
/// BSDFs since they may contain a diffuse component and a perfect specular
/// component. If the base path is currently working with a sample from a BSDF's
/// smooth component, we don't want to care about the specular component of the
/// BSDF right now - we want to deal with the smooth component.
///
/// For this reason, we vary the classification a little bit based on the
/// situation. This is perfectly valid, and should be done.
inline VertexType getVertexType(const BSDF *bsdf, Intersection &its,
                                const GradientPathTracerConfig &config,
                                unsigned int                    bsdfType) {
  // Return the lowest roughness value of the components of the vertex's BSDF.
  // If 'bsdfType' does not have a delta component, do not take perfect
  // speculars (zero roughness) into account in this.

  Float lowest_roughness = std::numeric_limits<Float>::infinity();

  bool found_smooth = false;
  bool found_dirac  = false;
  for (int i = 0, component_count = bsdf->getComponentCount();
       i < component_count; ++i) {
    Float component_roughness = bsdf->getRoughness(its, i);

    if (component_roughness == Float(0)) {
      found_dirac = true;
      if (!(bsdfType & BSDF::EDelta)) {
        // Skip Dirac components if a smooth component is requested.
        continue;
      }
    } else {
      found_smooth = true;
    }

    if (component_roughness < lowest_roughness) {
      lowest_roughness = component_roughness;
    }
  }

  // Roughness has to be zero also if there is a delta component but no smooth
  // components.
  if (!found_smooth && found_dirac && !(bsdfType & BSDF::EDelta)) {
    lowest_roughness = Float(0);
  }

  return getVertexTypeByRoughness(lowest_roughness, config);
}

inline VertexType getVertexType(RayState                       &ray,
                                const GradientPathTracerConfig &config,
                                unsigned int                    bsdfType) {
  const BSDF *bsdf = ray.rRec.its.getBSDF(ray.ray);
  return getVertexType(bsdf, ray.rRec.its, config, bsdfType);
}

/**
 * ------------------------- GPT Interface -------------------------
 */
class IGradientPathTracer {
public:
  IGradientPathTracer(const Scene *scene, const Sensor *sensor,
                      Sampler *sampler, GPTWorkResult *block,
                      const GradientPathTracerConfig *config)
      : m_scene(scene), m_sensor(sensor), m_sampler(sampler), m_block(block),
        m_config(config) {}

  virtual ~IGradientPathTracer() = default;

  /// Evaluates a sample at the given position.
  ///
  /// Outputs direct radiance to be added on top of the final image, the
  /// throughput to the central pixel, gradients to all neighbors, and
  /// throughput contribution to the neighboring pixels.
  virtual void evaluatePoint(RadianceQueryRecord &rRec,
                             const Point2        &samplePosition,
                             const Point2 &apertureSample, Float timeSample,
                             Float     differentialScaleFactor,
                             Spectrum &out_very_direct,
                             Spectrum &out_throughput, Spectrum *out_gradients,
                             Spectrum *out_neighborThroughputs) = 0;

  /// Constructs a sequence of base paths and shifts them into offset paths,
  /// evaluating their throughputs and differences.
  ///
  /// This is the core of the rendering algorithm.
  virtual void evaluate(RayState &main, RayState *shiftedRays,
                        int secondaryCount, Spectrum &out_veryDirect) = 0;

  /// Samples a direction according to the BSDF at the given ray position.
  inline BSDFSampleResult sampleBSDF(RayState &rayState) {
    Intersection        &its  = rayState.rRec.its;
    RadianceQueryRecord &rRec = rayState.rRec;
    RayDifferential     &ray  = rayState.ray;

    // Note: If the base path's BSDF evaluation uses random numbers, it would be
    // beneficial to use the same random numbers for the offset path's BSDF.
    //       This is not done currently.

    const BSDF *bsdf = its.getBSDF(ray);

    // Sample BSDF * cos(theta).
    BSDFSampleResult result = {BSDFSamplingRecord(its, rRec.sampler, ERadiance),
                               Spectrum(), (Float)0};

    Point2 sample = rRec.nextSample2D();
    result.weight = bsdf->sample(result.bRec, result.pdf, sample);

    // Variable result.pdf will be 0 if the BSDF sampler failed to produce a
    // valid direction.

    SAssert(result.pdf <= (Float)0 ||
            fabs(result.bRec.wo.length() - 1.0) < 0.00001);
    return result;
  }

protected:
  const Scene                    *m_scene;
  const Sensor                   *m_sensor;
  Sampler                        *m_sampler;
  GPTWorkResult                  *m_block;
  const GradientPathTracerConfig *m_config;
};

MTS_NAMESPACE_END
