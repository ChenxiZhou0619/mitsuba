#include <mitsuba/core/statistics.h>
#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

// TODO subsurface

class SpectralPathTracer : public SpectralMonteCarloIntegrator {
public:
  SpectralPathTracer(const Properties &props) : SpectralMonteCarloIntegrator(props) {}

  SpectralPathTracer(Stream *stream, InstanceManager *manager)
      : SpectralMonteCarloIntegrator(stream, manager) {}

  SampledSpectrum Li(const RayDifferential &r, SampledWavelengths &lambda,
                     RadianceQueryRecord &rRec) const {
    const Scene    *scene = rRec.scene;
    Intersection   &its   = rRec.its;
    RayDifferential ray(r);

    SampledSpectrum Li(.0f);
    SampledSpectrum beta(1.f);

    bool scattered = false;

    rRec.rayIntersect(ray);
    ray.mint = Epsilon;

    Float eta     = 1.f;
    int   bounces = 0;

    while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
      if (!its.isValid()) {
        /* If no intersection could be found, potentially return
           radiance from a environment luminaire if it exists */
        if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && (!m_hideEmitters || scattered))
          Li += beta * RGBIlluminantToSampledSpectrum(scene->evalEnvironment(ray), lambda);
        break;
      }

      const BSDF *bsdf = its.getBSDF(ray);

      /* Possibly include emitted radiance if requested */
      if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) &&
          (!m_hideEmitters || scattered))
        Li += beta * RGBIlluminantToSampledSpectrum(its.Le(-ray.d), lambda);

      if ((rRec.depth >= m_maxDepth && m_maxDepth > 0) ||
          (m_strictNormals && dot(ray.d, its.geoFrame.n) * Frame::cosTheta(its.wi) >= 0)) {

        /* Only continue if:
           1. The current path length is below the specifed maximum
           2. If 'strictNormals'=true, when the geometric and shading
              normals classify the incident direction to the same side */
        break;
      }

      /* ==================================================================== */
      /*                     Direct illumination sampling                     */
      /* ==================================================================== */

      /* Estimate the direct illumination if this is requested */
      DirectSamplingRecord dRec(its);

      if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
          (bsdf->getType() & BSDF::ESmooth)) {
        SampledSpectrum value = RGBIlluminantToSampledSpectrum(
            scene->sampleEmitterDirect(dRec, rRec.nextSample2D()), lambda);
        if (value) {
          const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

          /* Allocate a record for querying the BSDF */
          BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

          /* Evaluate BSDF * cos(theta) */
          SampledSpectrum bsdfVal = RGBAlbedoToSampledSpectrum(bsdf->eval(bRec), lambda);

          /* Prevent light leaks due to the use of shading normals */
          if (bsdfVal &&
              (!m_strictNormals || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

            /* Calculate prob. of having generated that direction
               using BSDF sampling */
            Float bsdfPdf =
                (emitter->isOnSurface() && dRec.measure == ESolidAngle) ? bsdf->pdf(bRec) : 0;

            /* Weight using the power heuristic */
            Float weight = miWeight(dRec.pdf, bsdfPdf);
            Li += beta * value * bsdfVal * weight;
          }
        }
      }

      /* ==================================================================== */
      /*                            BSDF sampling                             */
      /* ==================================================================== */

      /* Sample BSDF * cos(theta) */
      Float              bsdfPdf;
      BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);

      SampledSpectrum bsdfWeight =
          RGBAlbedoToSampledSpectrum(bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D()), lambda);
      if (!bsdfWeight) break;

      scattered |= bRec.sampledType != BSDF::ENull;

      /* Prevent light leaks due to the use of shading normals */
      const Vector wo        = its.toWorld(bRec.wo);
      Float        woDotGeoN = dot(its.geoFrame.n, wo);
      if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0) break;

      bool            hitEmitter = false;
      SampledSpectrum value;

      /* Trace a ray in this direction */
      ray = Ray(its.p, wo, ray.time);
      if (scene->rayIntersect(ray, its)) {
        /* Intersected something - check if it was a luminaire */
        if (its.isEmitter()) {
          value = RGBIlluminantToSampledSpectrum(its.Le(-ray.d), lambda);
          dRec.setQuery(ray, its);
          hitEmitter = true;
        }
      } else {
        /* Intersected nothing -- perhaps there is an environment map? */
        const Emitter *env = scene->getEnvironmentEmitter();

        if (env) {
          if (m_hideEmitters && !scattered) break;

          value = RGBIlluminantToSampledSpectrum(env->evalEnvironment(ray), lambda);
          if (!env->fillDirectSamplingRecord(dRec, ray)) break;
          hitEmitter = true;
        } else {
          break;
        }
      }

      /* Keep track of the throughput and relative
         refractive index along the path */
      beta *= bsdfWeight;
      eta *= bRec.eta;

      /* If a luminaire was hit, estimate the local illumination and
         weight using the power heuristic */
      if (hitEmitter && (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
        /* Compute the prob. of generating that direction using the
           implemented direct illumination sampling technique */
        const Float lumPdf =
            (!(bRec.sampledType & BSDF::EDelta)) ? scene->pdfEmitterDirect(dRec) : 0;
        Li += beta * value * miWeight(bsdfPdf, lumPdf);
      }

      /* ==================================================================== */
      /*                         Indirect illumination                        */
      /* ==================================================================== */

      /* Set the recursive query type. Stop if no surface was hit by the
         BSDF sample or if indirect illumination was not requested */
      if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) break;
      rRec.type = RadianceQueryRecord::ERadianceNoEmission;

      if (rRec.depth++ >= m_rrDepth) {
        /* Russian roulette: try to keep path weights equal to one,
           while accounting for the solid angle compression at refractive
           index boundaries. Stop with at least some probability to avoid
           getting stuck (e.g. due to total internal reflection) */

        Float q = std::min(beta.max() * eta * eta, (Float)0.95f);
        if (rRec.nextSample1D() >= q) break;
        beta /= q;
      }
    }

    return Li;
  }

protected:
  Float miWeight(Float pdfA, Float pdfB) const {
    return (pdfA * pdfA) / (pdfA * pdfA + pdfB * pdfB);
  }

  MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(SpectralPathTracer, false, SpectralMonteCarloIntegrator)
MTS_EXPORT_PLUGIN(SpectralPathTracer, "Spectral hero-wavelength path tracer")

MTS_NAMESPACE_END