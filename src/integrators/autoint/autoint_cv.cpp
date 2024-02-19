#include "iilf.h"
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/scene.h>
MTS_NAMESPACE_BEGIN

struct TrainingInfo {
  Point2i pixel_offset;

  Vector wi_local;
  Spectrum li;
};

class AutointControlVariate : public MonteCarloIntegrator {
public:
  AutointControlVariate(const Properties &props) : MonteCarloIntegrator(props) {
    m_only_indirect = props.getBoolean("only_indirect", true);
    m_hide_emitters = props.getBoolean("hide_emitters", true);
    m_training_samples = props.getInteger("training_samples", 8);
  }

  AutointControlVariate(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {}

  Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
    //! Unused
    return Spectrum(.5f);
  }

  Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec,
              TrainingInfo *info, std::shared_ptr<IILF> iilf) const {

    /* Some aliases and local variables */
    const Scene *scene = rRec.scene;
    Intersection &its = rRec.its;
    RayDifferential ray(r);
    Spectrum Li(0.0f);
    bool scattered = false;

    /* Perform the first ray intersection (or ignore if the
       intersection has already been provided). */
    rRec.rayIntersect(ray);
    ray.mint = Epsilon;

    Spectrum throughput(1.0f);
    Float eta = 1.0f;

    while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
      if (!its.isValid()) {
        /* If no intersection could be found, potentially return
           radiance from a environment luminaire if it exists */
        if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) &&
            (!m_hide_emitters || scattered))
          Li += throughput * scene->evalEnvironment(ray);
        break;
      }

      const BSDF *bsdf = its.getBSDF(ray);

      /* Possibly include emitted radiance if requested */
      if (its.isEmitter() &&
          (rRec.type & RadianceQueryRecord::EEmittedRadiance) &&
          (!m_hide_emitters || scattered))
        Li += throughput * its.Le(-ray.d);

      /* Include radiance from a subsurface scattering model if requested */
      if (its.hasSubsurface() &&
          (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
        Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

      if ((rRec.depth >= m_maxDepth && m_maxDepth > 0) ||
          (m_strictNormals &&
           dot(ray.d, its.geoFrame.n) * Frame::cosTheta(its.wi) >= 0)) {

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
          (bsdf->getType() & BSDF::ESmooth) &&
          (!m_only_indirect || rRec.depth != 1)) {
        Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
        if (!value.isZero()) {
          const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

          /* Allocate a record for querying the BSDF */
          BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

          /* Evaluate BSDF * cos(theta) */
          const Spectrum bsdfVal = bsdf->eval(bRec);

          /* Prevent light leaks due to the use of shading normals */
          if (!bsdfVal.isZero() &&
              (!m_strictNormals ||
               dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

            /* Calculate prob. of having generated that direction
               using BSDF sampling */
            Float bsdfPdf =
                (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                    ? bsdf->pdf(bRec)
                    : 0;

            /* Weight using the power heuristic */
            Float weight = miWeight(dRec.pdf, bsdfPdf);
            Li += throughput * value * bsdfVal * weight;
          }
        }
      }

      /* ==================================================================== */
      /*                            BSDF sampling                             */
      /* ==================================================================== */

      /* Sample BSDF * cos(theta) */
      Float bsdfPdf;
      BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
      Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
      if (bsdfWeight.isZero())
        break;

      scattered |= bRec.sampledType != BSDF::ENull;

      /* Prevent light leaks due to the use of shading normals */
      const Vector wo = its.toWorld(bRec.wo);
      Float woDotGeoN = dot(its.geoFrame.n, wo);
      if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
        break;

      bool hitEmitter = false;
      Spectrum value;

      /* Trace a ray in this direction */
      ray = Ray(its.p, wo, ray.time);
      if (scene->rayIntersect(ray, its)) {
        /* Intersected something - check if it was a luminaire */
        if (its.isEmitter()) {
          value = its.Le(-ray.d);
          dRec.setQuery(ray, its);
          hitEmitter = true;
        }
      } else {
        /* Intersected nothing -- perhaps there is an environment map? */
        const Emitter *env = scene->getEnvironmentEmitter();

        if (env) {
          if (m_hide_emitters && !scattered)
            break;

          value = env->evalEnvironment(ray);
          if (!env->fillDirectSamplingRecord(dRec, ray))
            break;
          hitEmitter = true;
        } else {
          break;
        }
      }

      /* Keep track of the throughput and relative
         refractive index along the path */
      throughput *= bsdfWeight;
      eta *= bRec.eta;

      /* If a luminaire was hit, estimate the local illumination and
         weight using the power heuristic */
      if (hitEmitter &&
          (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) &&
          (!m_only_indirect || rRec.depth != 1)) {
        /* Compute the prob. of generating that direction using the
           implemented direct illumination sampling technique */
        const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta))
                                 ? scene->pdfEmitterDirect(dRec)
                                 : 0;
        Li += throughput * value * miWeight(bsdfPdf, lumPdf);
      }

      /* ==================================================================== */
      /*                         Indirect illumination                        */
      /* ==================================================================== */

      /* Set the recursive query type. Stop if no surface was hit by the
         BSDF sample or if indirect illumination was not requested */
      if (!its.isValid() ||
          !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
        break;
      rRec.type = RadianceQueryRecord::ERadianceNoEmission;

      if (rRec.depth++ >= m_rrDepth) {
        /* Russian roulette: try to keep path weights equal to one,
           while accounting for the solid angle compression at refractive
           index boundaries. Stop with at least some probability to avoid
           getting stuck (e.g. due to total internal reflection) */

        Float q = std::min(throughput.max() * eta * eta, (Float)0.95f);
        if (rRec.nextSample1D() >= q)
          break;
        throughput /= q;
      }
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
    bool needsTimeSample = sensor->needsTimeSample();

    RadianceQueryRecord rRec(scene, sampler);
    Point2 apertureSample(0.5f);
    Float timeSample = 0.5f;
    RayDifferential sensorRay;

    block->clear();

    uint32_t queryType = RadianceQueryRecord::ESensorRay;

    if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we
                                           don't have to */
      queryType &= ~RadianceQueryRecord::EOpacity;

    for (size_t i = 0; i < points.size(); ++i) {
      Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
      int linear_offset = offset[0] * 288 + offset[1]; // TODO

      if (stop)
        break;

      sampler->generate(offset);

      std::shared_ptr<IILF> iilf = nullptr;
      if (!m_is_training /* real rendering */) {
        // fit the autoint with samples stored in m_datas

        iilf = std::make_shared<IILF>();
        iilf->fit(m_datas[linear_offset]);
      }

      for (size_t j = 0; j < sampler->getSampleCount(); j++) {
        rRec.newQuery(queryType, sensor->getMedium());
        Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

        if (needsApertureSample)
          apertureSample = rRec.nextSample2D();
        if (needsTimeSample)
          timeSample = rRec.nextSample1D();

        Spectrum spec = sensor->sampleRayDifferential(
            sensorRay, samplePos, apertureSample, timeSample);

        sensorRay.scaleDifferential(diffScaleFactor);

        TrainingInfo info{offset, Vector(.0f), Spectrum(.0f)};

        spec *= Li(sensorRay, rRec, &info, iilf);

        if (m_is_training) {
          m_datas[linear_offset].emplace_back(
              std::pair{info.wi_local, info.li});
        }

        block->put(samplePos, spec, rRec.alpha);
        sampler->advance();
      }
    }
  }

  MTS_DECLARE_CLASS()

private:
  bool m_only_indirect;
  bool m_hide_emitters;
  int m_training_samples;

  bool m_is_training;

  static constexpr int m_screen_size = 512 * 288;

  mutable std::vector<std::pair<Vector, Spectrum>> m_datas[m_screen_size];
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
  ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
  ref<Scene> train_scene = new Scene(scene);

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