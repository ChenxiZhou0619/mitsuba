// clang-format off
#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>
#include "cache/field.h"
#include <mitsuba/core/plugin.h>
// clang-format on

MTS_NAMESPACE_BEGIN

// KD-tree + spherical harmonic radiance cache
// Only surface

class RadianceCachePathTracer : public MonteCarloIntegrator {
public:
  RadianceCachePathTracer(const Properties &props)
      : MonteCarloIntegrator(props) {
    //
  }

  RadianceCachePathTracer(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {
    Log(EError, "RadianceCache PT serialization is not support");
  }

  void serialize(Stream *stream, InstanceManager *manager) const override {
    MonteCarloIntegrator::serialize(stream, manager);
    Log(EError, "RadianceCache PT serialization is not support");
  }

  virtual bool preprocess(const Scene *scene, RenderQueue *queue,
                          const RenderJob *job, int sceneResID, int sensorResID,
                          int samplerResID) override;

  virtual Spectrum Li(const RayDifferential &r,
                      RadianceQueryRecord   &rRec) const override {

    const Scene    *scene = rRec.scene;
    Intersection   &its   = rRec.its;
    RayDifferential ray(r);

    Spectrum L(.0f), beta(1.f);
    bool     specularBounce = false;

    Float  prevScatterPDF = .0f;
    Point  prevP;
    Normal prevN;

    while (true) {
      if (beta.isZero()) break;

      bool foundIntersection = scene->rayIntersect(ray, its);

      if (!foundIntersection) {
        Spectrum Le = scene->evalEnvironment(ray);
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
        Spectrum Le = emitter->eval(its, -ray.d);
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
        L += beta * SampleNEE(scene, dRec, -ray.d, rRec.sampler, bsdf, &its);

      BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
      Float              bsdfPDF;
      Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPDF, rRec.nextSample2D());
      if (bsdfWeight.isZero()) break;

      const Vector wo        = its.toWorld(bRec.wo);
      Float        woDotGeoN = dot(its.geoFrame.n, wo);
      if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals) break;

      beta *= bsdfWeight;
      ray = Ray(its.p + Epsilon * wo, wo, Epsilon, INFINITY, ray.time);
      specularBounce = (bRec.sampledType &
                        (BSDF::EDeltaReflection | BSDF::EDeltaTransmission));

      prevScatterPDF = bsdfPDF;
      prevP          = its.p;
      prevN          = its.shFrame.n;

      if (rRec.depth >= m_rrDepth) {
        Float q = std::min(beta.max(), (Float)0.95f);
        if (rRec.nextSample1D() >= q) break;
        beta /= q;
      }
    }

    return L;
  }

protected:
  void train(Scene *scene, RenderQueue *queue, const RenderJob *job,
             int sensorResID);

  Spectrum SampleNEE(const Scene *scene, DirectSamplingRecord &dRec,
                     const Vector &wi, Sampler *sampler, const BSDF *bsdf,
                     const Intersection *its) const {
    Spectrum Le(.0f);
    Spectrum scatterVal(.0f);
    Float    scatterPDF = .0f;
    Point2   sample     = sampler->next2D();

    Le = scene->sampleEmitterDirect(dRec, sample, true);
    if (Le.isZero()) return Spectrum(.0f);

    const Emitter     *emitter = static_cast<const Emitter *>(dRec.object);
    BSDFSamplingRecord bRec(*its, its->toLocal(dRec.d), ERadiance);
    scatterVal = bsdf->eval(bRec);
    scatterPDF = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                     ? bsdf->pdf(bRec)
                     : .0f;
    Float misw = MISWeight(dRec.pdf, scatterPDF);
    return scatterVal * Le * misw;
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

  Float MISWeight(Float a, Float b) const { return a / (a + b); }

  MTS_DECLARE_CLASS()

private:
  // gather the radiance data at each
  uint32_t                        m_num_iterations = 3;
  bool                            m_training       = false;
  std::unique_ptr<RadianceKDTree> m_radianceCache;
  std::vector<RadianceRecord>     m_data;
};

bool RadianceCachePathTracer::preprocess(const Scene *scene, RenderQueue *queue,
                                         const RenderJob *job, int sceneResID,
                                         int sensorResID, int samplerResID) {

  m_radianceCache->setAABB(scene->getAABB());
  ref<Scheduler> sched = Scheduler::getInstance();
  train(static_cast<Scene *>(sched->getResource(sceneResID)), queue, job,
        sensorResID);
  return true;
}

void RadianceCachePathTracer::train(Scene *scene, RenderQueue *queue,
                                    const RenderJob *job, int sensorResID) {
  m_training = true;

  ref<Scheduler> sched = Scheduler::getInstance();
  ref<Sensor> sensor   = static_cast<Sensor *>(sched->getResource(sensorResID));
  ref<Scene>  train_scene = new Scene(scene);

  const int trainSceneResID = sched->registerResource(train_scene);

  for (int i = 0; i < m_num_iterations; ++i) {
    Properties training_sampler_props = scene->getSampler()->getProperties();

    training_sampler_props.removeProperty("sampleCount");
    training_sampler_props.setSize("sampleCount", std::pow(2, i));

    ref<Sampler> train_sampler =
        static_cast<Sampler *>(PluginManager::getInstance()->createObject(
            MTS_CLASS(Sampler), training_sampler_props));
    train_sampler->configure();
    train_scene->setSampler(train_sampler);

    //* Create a sampler for each thread
    std::vector<SerializableObject *> samplers(sched->getCoreCount());
    for (size_t i = 0; i < sched->getCoreCount(); ++i) {
      ref<Sampler> cloned_sampler = train_sampler->clone();
      cloned_sampler->incRef();
      samplers[i] = cloned_sampler.get();
    }
    int trainingSamplerResID = sched->registerMultiResource(samplers);
    for (size_t i = 0; i < sched->getCoreCount(); ++i)
      samplers[i]->decRef();

    SamplingIntegrator::render(train_scene, queue, job, trainSceneResID,
                               sensorResID, trainingSamplerResID);
  }
}

MTS_IMPLEMENT_CLASS(RadianceCachePathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(RadianceCachePathTracer, "Path tracer with radiance cache")

MTS_NAMESPACE_END