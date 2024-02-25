#include "func_2d.h"
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/scene.h>
MTS_NAMESPACE_BEGIN

struct TrainingInfo {
  Point2i pixel_offset;

  Vector   wo_local[2];
  Spectrum li[2];
};

class AutointControlVariate : public MonteCarloIntegrator {
public:
  AutointControlVariate(const Properties &props) : MonteCarloIntegrator(props) {
    m_training_samples = props.getInteger("training_samples", 8);
    m_screen_size      = props.getInteger("screen_size", 921600);
  }

  AutointControlVariate(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {}

  Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
    //! Unused
    return Spectrum(.5f);
  }

  Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec,
              TrainingInfo *info, std::shared_ptr<Func_2d> func_2d) const {
    const Scene    *scene = rRec.scene;
    Intersection   &its   = rRec.its;
    RayDifferential ray(r);
    Spectrum        Li(.0f);
    bool            scattered = false;

    rRec.rayIntersect(ray);
    ray.mint = Epsilon;

    Spectrum throughput(1.f);
    Float    eta = 1.0f;

    if (!its.isValid()) {
      if (!m_hideEmitters) Li += throughput * scene->evalEnvironment(ray);

      return Li;
    }

    const BSDF *bsdf = its.getBSDF(ray);

    if (its.isEmitter() && !m_hideEmitters) Li += throughput * its.Le(-ray.d);

    if ((m_strictNormals &&
         dot(ray.d, its.geoFrame.n) * Frame::cosTheta(its.wi) >= 0))
      return Li;

    //* Direct illumination

    DirectSamplingRecord dRec(its);

    //* NEE
    {
      Vector   wo_local;
      Spectrum li(.0f);
      if (bsdf->getType() & BSDF::ESmooth) {
        Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
        if (!value.isZero()) {

          const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

          BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

          const Spectrum bsdf_val = bsdf->eval(bRec);

          wo_local = bRec.wo;
          if (!bsdf_val.isZero() &&
              dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0) {
            Float bsdf_pdf =
                (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                    ? bsdf->pdf(bRec)
                    : .0f;
            Float weight = miWeight(dRec.pdf, bsdf_pdf);
            li = bsdf_val * value * dRec.pdf * Frame::sinTheta(wo_local);

            if (m_is_training)
              Li += throughput * value * bsdf_val * weight;
            else {
              auto [analytic_int, cv] = func_2d->control_variate(bRec.wo);

              Spectrum fx = dRec.pdf * value * bsdf_val;
              Spectrum gx = cv / Frame::sinTheta(bRec.wo);

              Spectrum estimator = analytic_int + (fx - gx) / dRec.pdf;

              Li += throughput * weight * estimator;
            }
          }
        }
      }

      info->wo_local[0] = wo_local;
      info->li[0]       = li;
    }

    { //* BSDF Sampling
      Vector   wo_local;
      Spectrum li(.0f);

      Float              bsdf_pdf;
      BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
      Spectrum bsdf_weight = bsdf->sample(bRec, bsdf_pdf, rRec.nextSample2D());
      if (bsdf_weight.isZero()) return Li;

      const Vector wo        = its.toWorld(bRec.wo);
      Float        woDotGeoN = dot(its.geoFrame.n, wo);

      wo_local = bRec.wo;

      if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
        return Li;

      bool     hit_emitter = false;
      Spectrum value;

      ray = Ray(its.p, wo, ray.time);
      if (scene->rayIntersect(ray, its)) {
        if (its.isEmitter()) {
          value = its.Le(-ray.d);
          dRec.setQuery(ray, its);
          hit_emitter = true;
        }
      } else {
        const Emitter *env = scene->getEnvironmentEmitter();

        if (env) {
          if (m_hideEmitters) return Li;

          value = env->evalEnvironment(ray);
          if (!env->fillDirectSamplingRecord(dRec, ray)) return Li;
          hit_emitter = true;
        } else {
          return Li;
        }
      }

      eta *= bRec.eta;

      if (hit_emitter) {
        const Float lum_pdf = (!(bRec.sampledType & BSDF::EDelta))
                                  ? scene->pdfEmitterDirect(dRec)
                                  : 0;

        li         = bsdf_weight * bsdf_pdf * Frame::sinTheta(bRec.wo) * value;
        Float misw = miWeight(bsdf_pdf, lum_pdf);

        if (m_is_training)
          Li += throughput * bsdf_weight * value * misw;
        else {
          auto [analytic_int, cv] = func_2d->control_variate(bRec.wo);

          Spectrum fx = bsdf_weight * bsdf_pdf * value;
          Spectrum gx = cv / Frame::sinTheta(bRec.wo);

          Spectrum estimator = analytic_int + (fx - gx) / bsdf_pdf;

          Li += throughput * misw * estimator;
        }
      }

      info->wo_local[1] = wo_local;
      info->li[1]       = li;
    }
    return Li;
  }

  inline Float miWeight(Float pdfA, Float pdfB) const {
    pdfA *= pdfA;
    pdfB *= pdfB;
    return pdfA / (pdfA + pdfB);
  }

  void serialize(Stream *stream, InstanceManager *manager) const {
    MonteCarloIntegrator::serialize(stream, manager);
    std::cout << "AutointControlVariate::serialize not implemented!\n";
  }

  std::string toString() const {
    return "AutointControlVariate screen-spaced path tracer\n";
  }

  bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
                  int sceneResID, int sensorResID, int samplerResID);

  void train(Scene *scene, RenderQueue *queue, const RenderJob *job,
             int sensorResID);

  void renderBlock(const Scene *scene, const Sensor *sensor, Sampler *sampler,
                   ImageBlock *block, const bool &stop,
                   const std::vector<TPoint2<uint8_t>> &points) const {

    Float diffScaleFactor = 1.0f / std::sqrt((Float)sampler->getSampleCount());

    bool needsApertureSample = sensor->needsApertureSample();
    bool needsTimeSample     = sensor->needsTimeSample();

    RadianceQueryRecord rRec(scene, sampler);
    Point2              apertureSample(0.5f);
    Float               timeSample = 0.5f;
    RayDifferential     sensorRay;

    block->clear();

    uint32_t queryType = RadianceQueryRecord::ESensorRay;

    if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we
                                           don't have to */
      queryType &= ~RadianceQueryRecord::EOpacity;

    for (size_t i = 0; i < points.size(); ++i) {
      Point2i offset        = Point2i(points[i]) + Vector2i(block->getOffset());
      int     linear_offset = offset[0] * 720 + offset[1]; // TODO

      if (stop) break;

      sampler->generate(offset);

      std::shared_ptr<Func_2d> func_2d = nullptr;
      if (!m_is_training /* real rendering */) {
        // fit the autoint with samples stored in m_datas

        func_2d = std::make_shared<Func_2d>();
        func_2d->fit(m_datas[linear_offset]);
      }

      for (size_t j = 0; j < sampler->getSampleCount(); j++) {
        rRec.newQuery(queryType, sensor->getMedium());
        Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

        if (needsApertureSample) apertureSample = rRec.nextSample2D();
        if (needsTimeSample) timeSample = rRec.nextSample1D();

        Spectrum spec = sensor->sampleRayDifferential(
            sensorRay, samplePos, apertureSample, timeSample);

        sensorRay.scaleDifferential(diffScaleFactor);

        TrainingInfo info{offset};

        spec *= Li(sensorRay, rRec, &info, func_2d);

        if (m_is_training) {
          m_datas[linear_offset].emplace_back(
              std::pair{info.wo_local[0], info.li[0]});
          m_datas[linear_offset].emplace_back(
              std::pair{info.wo_local[1], info.li[1]});
        }

        block->put(samplePos, spec, rRec.alpha);
        sampler->advance();
      }
    }
  }

  MTS_DECLARE_CLASS()

private:
  int  m_training_samples, m_screen_size;
  bool m_is_training;

  // static constexpr int m_screen_size = 512 * 288;

  mutable std::vector<std::pair<Vector, Spectrum>> m_datas[921600];
};

bool AutointControlVariate::preprocess(const Scene *scene, RenderQueue *queue,
                                       const RenderJob *job, int sceneResID,
                                       int sensorResID, int samplerResID) {
  ref<Scheduler> sched = Scheduler::getInstance();
  train(static_cast<Scene *>(sched->getResource(sceneResID)), queue, job,
        sensorResID);
  return true;
}

void AutointControlVariate::train(Scene *scene, RenderQueue *queue,
                                  const RenderJob *job, int sensorResID) {
  m_is_training = true;

  ref<Scheduler> sched = Scheduler::getInstance();
  ref<Sensor> sensor   = static_cast<Sensor *>(sched->getResource(sensorResID));
  ref<Scene>  train_scene = new Scene(scene);

  const int trainSceneResID = sched->registerResource(train_scene);

  int num_iterations = 1; // TODO
  for (int i = 0; i < num_iterations; ++i) {
    Properties training_sampler_props = scene->getSampler()->getProperties();

    training_sampler_props.removeProperty("sampleCount");
    training_sampler_props.setSize("sampleCount", m_training_samples);

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
    sched->unregisterResource(trainingSamplerResID);
    sensor->getFilm()->clear();
  }

  m_is_training = false;
}

MTS_IMPLEMENT_CLASS_S(AutointControlVariate, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(AutointControlVariate, "Autoint control variate pt");
MTS_NAMESPACE_END