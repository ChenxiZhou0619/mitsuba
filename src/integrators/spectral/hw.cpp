// clang-format off
#include "spectralintegrator.h"
#include <mitsuba/render/scene.h>

// clang-format on

MTS_NAMESPACE_BEGIN

using namespace spectral;

class HeroWavelengthPathTracer : public SpectralMonteCarloIntegrator {
public:
  HeroWavelengthPathTracer(const Properties &props)
      : SpectralMonteCarloIntegrator(props) {
    //
  }

  HeroWavelengthPathTracer(Stream *stream, InstanceManager *manager)
      : SpectralMonteCarloIntegrator(stream, manager) {
    Log(EError, "SpectralIntegrator PT serialization is not support");
  }

  void serialize(Stream *stream, InstanceManager *manager) const override {
    SpectralMonteCarloIntegrator::serialize(stream, manager);
    Log(EError, "SpectralIntegrator PT serialization is not support");
  }

  virtual SampledSpectrum Li(const RayDifferential &r,
                             RadianceQueryRecord   &rRec,
                             SampledWavelengths    &lambdas) const override {
    const Scene    *scene = rRec.scene;
    Intersection   &its   = rRec.its;
    RayDifferential ray(r);

    SampledSpectrum L(.0f), beta(1.f);
    bool            specularBounce = false;

    Float  prevScatterPDF = .0f;
    Point  prevP;
    Normal prevN;

    while (true) {
      if (!beta) break;

      bool foundIntersection = scene->rayIntersect(ray, its);

      if (!foundIntersection) {
        SampledSpectrum Le = RGBToSampledSpectrum(scene->evalEnvironment(ray),
                                                  lambdas, EIlluminantRGB);
        if (rRec.depth == 1 || specularBounce) {
          // First hit or pure specular bounce
          L += beta * Le;
        } else {
          // Apply MIS
          Float neePDF = PDF_Nee(scene, prevP, prevN, ray);
          Float misw   = MISWeight(prevScatterPDF, neePDF);
          L += beta * Le * misw;
        }
        break;
      }

      if (const Emitter *emitter = its.shape->getEmitter(); emitter) {
        SampledSpectrum Le = RGBToSampledSpectrum(emitter->eval(its, -ray.d),
                                                  lambdas, EIlluminantRGB);
        if (rRec.depth == 1 || specularBounce) {
          L += beta * Le;
        } else {
          // Apply MIS
          Float neePDF = PDF_Nee(scene, prevP, prevN, ray, &its);
          Float misw   = MISWeight(prevScatterPDF, neePDF);
          L += beta * Le * misw;
        }
      }

      const BSDF *bsdf = its.getBSDF(ray);
      if (rRec.depth++ >= m_maxDepth) break;
      // Sample NEE
      DirectSamplingRecord dRec(its);
      if (bsdf->getType() & BSDF::ESmooth)
        L += beta *
             SampleNEE(scene, dRec, -ray.d, rRec.sampler, bsdf, &its, lambdas);

      BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
      Float              bsdfPDF;
      Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPDF, rRec.nextSample2D());
      if (bsdfWeight.isZero()) break;

      const Vector wo        = its.toWorld(bRec.wo);
      Float        woDotGeoN = dot(its.geoFrame.n, wo);
      if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals) break;

      beta *= RGBToSampledSpectrum(bsdfWeight, lambdas, EReflectanceRGB);
      ray = Ray(its.p + Epsilon * wo, wo, Epsilon, INFINITY, ray.time);
      specularBounce = (bRec.sampledType &
                        (BSDF::EDeltaReflection | BSDF::EDeltaTransmission));

      prevScatterPDF = bsdfPDF;
      prevP          = its.p;
      prevN          = its.shFrame.n;

      // FIXME: disable rr cuz rr should be after the uni contribution
    }
    return L;
  }

protected:
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

  Float MISWeight(Float a, Float b) const { return a / (a + b); }

  SampledSpectrum SampleNEE(const Scene *scene, DirectSamplingRecord &dRec,
                            const Vector &wi, Sampler *sampler,
                            const BSDF *bsdf, const Intersection *its,
                            const SampledWavelengths &lambdas) const {
    Spectrum rgbLe(.0f);
    Spectrum scatterVal(.0f);
    Float    scatterPDF = .0f;
    Point2   sample     = sampler->next2D();

    rgbLe = scene->sampleEmitterDirect(dRec, sample, true);
    if (rgbLe.isZero()) return SampledSpectrum(.0f);

    const Emitter     *emitter = static_cast<const Emitter *>(dRec.object);
    BSDFSamplingRecord bRec(*its, its->toLocal(dRec.d), ERadiance);
    scatterVal = bsdf->eval(bRec);
    scatterPDF = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                     ? bsdf->pdf(bRec)
                     : .0f;
    Float misw = MISWeight(dRec.pdf, scatterPDF);
    return RGBToSampledSpectrum(scatterVal, lambdas, EReflectanceRGB) *
           RGBToSampledSpectrum(rgbLe, lambdas, EIlluminantRGB) * misw;
  }
};

MTS_NAMESPACE_END
