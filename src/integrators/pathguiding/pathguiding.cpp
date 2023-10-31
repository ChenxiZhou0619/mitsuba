// clang-format off

#include <mitsuba/core/statistics.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

#include <openpgl/cpp/OpenPGL.h>
// clang-format on
MTS_NAMESPACE_BEGIN

using PGL_Device              = openpgl::cpp::Device;
using PGL_Field               = openpgl::cpp::Field;
using PGL_SurfaceDistribution = openpgl::cpp::SurfaceSamplingDistribution;
using PGL_SampleStorage       = openpgl::cpp::SampleStorage;
using PGL_SampleData          = openpgl::cpp::SampleData;
using PGL_PathSegmentStorge   = openpgl::cpp::PathSegmentStorage;
using PGL_PathSegment         = openpgl::cpp::PathSegment;

class PathGuidingTracer : public MonteCarloIntegrator {
public:
  // TODO
  PathGuidingTracer(const Properties &props) : MonteCarloIntegrator(props) {

    //* Default SD structure KDTree + VMM
    m_device = new PGL_Device(PGL_DEVICE_TYPE_CPU_4);
    PGLFieldArguments fieldArgs;
    pglFieldArgumentsSetDefaults(fieldArgs, PGL_SPATIAL_STRUCTURE_KDTREE,
                                 PGL_DIRECTIONAL_DISTRIBUTION_QUADTREE);

    m_field = new PGL_Field(m_device, fieldArgs);
  }

  PathGuidingTracer(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {
    Log(EError, "PathGuidingTracer serialization is not support\n");
  }

  void serializa(Stream *stream, InstanceManager *manager) {
    MonteCarloIntegrator::serialize(stream, manager);
    Log(EError, "PathGuidingTracer serialization is not support\n");
  }

  virtual ~PathGuidingTracer() {
    delete m_field;
    delete m_device;
  }

  //* Perform spatio-directional distribution training
  virtual bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
                          int sceneResID, int sensorResID, int samplerResID) override;

  void train(Scene *scene, RenderQueue *queue, const RenderJob *job, int sensorResID);

  //* Perform pixel estimating by tracing a ray
  virtual Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const override {
    const Scene    *scene = rRec.scene;
    Intersection   &its   = rRec.its;
    RayDifferential ray(r);

    Spectrum Li(.0f);
    Spectrum beta(1.f);
    Float    eta = 1.f;

    bool scattered = false;

    rRec.rayIntersect(ray);
    ray.mint = Epsilon;

    //!
    PGL_PathSegmentStorge path_storge;
    path_storge.Reserve(m_maxDepth);

    while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
      if (!its.isValid()) {
        /* If no intersection could be found, potentially return
           radiance from a environment luminaire if it exists */
        if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && (!m_hideEmitters || scattered))
          Li += beta * scene->evalEnvironment(ray);
        break;
      }

      const BSDF *bsdf = its.getBSDF(ray);

      /* Possibly include emitted radiance if requested */
      if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) &&
          (!m_hideEmitters || scattered))
        Li += beta * its.Le(-ray.d);

      /* Include radiance from a subsurface scattering model if requested */
      if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
        Li += beta * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

      if ((rRec.depth >= m_maxDepth && m_maxDepth > 0) ||
          (m_strictNormals && dot(ray.d, its.geoFrame.n) * Frame::cosTheta(its.wi) >= 0)) {

        /* Only continue if:
           1. The current path length is below the specifed maximum
           2. If 'strictNormals'=true, when the geometric and shading
              normals classify the incident direction to the same side */
        break;
      }

      PGL_PathSegment path_segment;

      path_segment.position     = {its.p.x, its.p.y, its.p.z};
      path_segment.directionOut = {-ray.d.x, -ray.d.y, -ray.d.z};
      path_segment.normal       = {its.geoFrame.n.x, its.geoFrame.n.y, its.geoFrame.n.z};

      /* ==================================================================== */
      /*                     Direct illumination sampling                     */
      /* ==================================================================== */

      PGL_SurfaceDistribution distr(m_field);
      if (!m_training) {
        Float f           = (Float)rand() / RAND_MAX;
        bool  distr_valid = distr.Init(m_field, {its.p.x, its.p.y, its.p.z}, f);
      }

      /* Estimate the direct illumination if this is requested */
      DirectSamplingRecord dRec(its);

      if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
          (bsdf->getType() & BSDF::ESmooth)) {
        Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
        if (!value.isZero()) {
          const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

          /* Allocate a record for querying the BSDF */
          BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

          /* Evaluate BSDF * cos(theta) */
          const Spectrum bsdfVal = bsdf->eval(bRec);

          /* Prevent light leaks due to the use of shading normals */
          if (!bsdfVal.isZero() &&
              (!m_strictNormals || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

            /* Calculate prob. of having generated that direction
               using BSDF sampling */

            Float bsdfPdf =
                (emitter->isOnSurface() && dRec.measure == ESolidAngle) ? bsdf->pdf(bRec) : 0;

            if (!m_training) bsdfPdf = distr.PDF({dRec.d[0], dRec.d[1], dRec.d[2]});

            /* Weight using the power heuristic */
            Float    weight = miWeight(dRec.pdf, bsdfPdf);
            Spectrum Ld     = value * bsdfVal;
            Li += weight * beta * Ld;

            path_segment.directContribution = {Ld[0], Ld[1], Ld[2]};
            path_segment.miWeight           = weight;
          }
        }
      }

      /* ==================================================================== */
      /*                            BSDF sampling                             */
      /* ==================================================================== */

      /* Sample BSDF * cos(theta) */
      Float              bsdfPdf;
      BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
      Spectrum           bsdfWeight(.0f);
      if (m_training)
        bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
      else {
        auto [u1, u2] = rRec.nextSample2D();
        pgl_vec3f dir;
        bsdfPdf    = distr.SamplePDF({u1, u2}, dir);
        bRec.wo    = {its.toLocal({dir.x, dir.y, dir.z})};
        bsdfWeight = bsdf->eval(bRec) / bsdfPdf;
      }

      if (bsdfWeight.isZero()) break;

      scattered |= bRec.sampledType != BSDF::ENull;

      /* Prevent light leaks due to the use of shading normals */
      const Vector wo        = its.toWorld(bRec.wo);
      Float        woDotGeoN = dot(its.geoFrame.n, wo);
      if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0) break;

      bool     hitEmitter = false;
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
          if (m_hideEmitters && !scattered) break;

          value = env->evalEnvironment(ray);
          if (!env->fillDirectSamplingRecord(dRec, ray)) break;
          hitEmitter = true;
        } else {
          break;
        }
      }

      /* Keep track of the beta and relative
         refractive index along the path */
      beta *= bsdfWeight;
      eta *= bRec.eta;

      path_segment.scatteringWeight = {bsdfWeight[0], bsdfWeight[1], bsdfWeight[2]};
      path_segment.pdfDirectionIn   = bsdfPdf;
      path_segment.directionIn      = {wo[0], wo[1], wo[2]};

      /* If a luminaire was hit, estimate the local illumination and
         weight using the power heuristic */
      if (hitEmitter && (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
        /* Compute the prob. of generating that direction using the
           implemented direct illumination sampling technique */
        const Float lumPdf =
            (!(bRec.sampledType & BSDF::EDelta)) ? scene->pdfEmitterDirect(dRec) : 0;
        Spectrum Ld = value * miWeight(bsdfPdf, lumPdf);
        Li += beta * Ld;
        path_segment.scatteredContribution = {Ld[0], Ld[1], Ld[2]};
      }

      if (m_training) {
        path_storge.AddSegment(path_segment);
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

    if (m_training) {
      size_t               n_samples = path_storge.PrepareSamples(true, true);
      const PGLSampleData *data      = path_storge.GetSamples(n_samples);
      m_storge.AddSamples(data, n_samples);
    }

    return Li;
  }

  inline Float miWeight(Float pdfA, Float pdfB) const {
    pdfA *= pdfA;
    pdfB *= pdfB;
    return pdfA / (pdfA + pdfB);
  }

private:
  //* Arguments

  PGL_Device               *m_device = nullptr;
  PGL_Field                *m_field  = nullptr;
  mutable PGL_SampleStorage m_storge;

  bool m_training;

  int m_n_iterations;
  int m_n_samples_per_itr;

  MTS_DECLARE_CLASS()
};

bool PathGuidingTracer::preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
                                   int sceneResID, int sensorResID, int samplerResID) {
  // set filed bound
  auto [x_min, y_min, z_min] = scene->getKDTree()->getAABB().min;
  auto [x_max, y_max, z_max] = scene->getKDTree()->getAABB().max;
  m_field->SetSceneBounds({{x_min, y_min, z_min}, {x_max, y_max, z_max}});

  ref<Scheduler> sched = Scheduler::getInstance();
  train(static_cast<Scene *>(sched->getResource(sceneResID)), queue, job, sensorResID);

  return true;
}

void PathGuidingTracer::train(Scene *scene, RenderQueue *queue, const RenderJob *job,
                              int sensorResID) {
  m_training = true;

  ref<Scheduler> sched       = Scheduler::getInstance();
  ref<Sensor>    sensor      = static_cast<Sensor *>(sched->getResource(sensorResID));
  ref<Scene>     train_scene = new Scene(scene);

  const int  trainSceneResID        = sched->registerResource(train_scene);
  Properties training_sampler_props = scene->getSampler()->getProperties();

  training_sampler_props.removeProperty("sampleCount");
  training_sampler_props.setSize("sampleCount", 16); // TODO

  ref<Sampler> train_sampler = static_cast<Sampler *>(
      PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), training_sampler_props));
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

  SamplingIntegrator::render(train_scene, queue, job, trainSceneResID, sensorResID,
                             trainingSamplerResID);

  m_field->Update(m_storge);
  m_field->Validate();

  Log(EInfo, "Training end");

  m_training = false;
  sched->unregisterResource(trainingSamplerResID);
}

MTS_IMPLEMENT_CLASS(PathGuidingTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(PathGuidingTracer, "Path guiding path tracing")

MTS_NAMESPACE_END
