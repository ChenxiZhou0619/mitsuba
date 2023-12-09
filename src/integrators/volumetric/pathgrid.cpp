// clang-format off

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

#include "pathgrid/pathgrid.h"
#include "pathgrid/path.h"

/**
 * A volumetric path-tracer using subpath reuse for multiple scattering
 * 
 * 
 * ! m_maxDepth + 1 == the maximum vertices number (include camera and lightsource)
 */

// clang-format on
MTS_NAMESPACE_BEGIN

class PathGrid : public MonteCarloIntegrator {
public:
  PathGrid(const Properties &props) : MonteCarloIntegrator(props) {
    m_training_spp = props.getInteger("training_spp", 16);
    m_training     = false;

    constexpr int cache_size = 1 << 20;

    m_pathGrid = std::make_unique<pathgrid::PathGrid>(cache_size);
    m_storage  = std::make_unique<pathgrid::PathStorage>(cache_size);
  }

  PathGrid(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {
    Log(EError, "PathGrid serialization is not support\n");
  }

  void serialize(Stream *stream, InstanceManager *manager) const override {
    MonteCarloIntegrator::serialize(stream, manager);
    Log(EError, "Path serialization is not support\n");
  }

  virtual Spectrum Li(const RayDifferential &r,
                      RadianceQueryRecord   &rRec) const override {
    using RNG = RNG_CX;

    // Uniformly sample a hero wave length
    const uint32_t heroChannel = SampleHeroChannel(rRec.nextSample1D());

    const Scene    *scene = rRec.scene;
    Intersection   &its   = rRec.its;
    RayDifferential ray(r);

    Spectrum L(.0f), beta(1.f), r_u(1.f), r_l(1.f);

    Float  prev_scatter_pdf = .0f;
    Point  prev_p;
    Normal prev_n;
    bool   specular_bounce = true;

    pathgrid::PathInfo path(m_maxDepth);
    path.addVertex(Point(.0f), Vector(.0f)); ///< dummy vertex

    while (true) {
      if (beta.isZero()) break;

      ray.heroChannel = heroChannel;
      scene->rayIntersect(ray, its);

      if (const Medium *medium = rRec.medium; medium) {
        bool scattered  = false;
        bool terminated = false;

        Float    t_max = its.isValid() ? its.t : INFINITY;
        uint64_t seed  = ((uint64_t)rand() << 32) + (uint64_t)rand();
        RNG      rng(seed);

        Spectrum T_maj = SampleTrMajorant(
            medium, ray, t_max, rng, ETrackingType::ENaiveDeltaTracking,
            [&](MajorantSamplingRecord majRec) {
              // Sample a type of collision
              Spectrum sigma_a = majRec.sigma_a;
              Spectrum sigma_s = majRec.sigma_s;
              Spectrum sigma_n = majRec.sigma_n;
              Spectrum tr_maj  = majRec.tr_majorant;

              ECollisionType scatter_type = SampleCollision(
                  sigma_a, sigma_s, sigma_n, rng.randomFloat(), heroChannel);

              if (scatter_type & ECollisionType::Absorb) {
                terminated = true;
                return false;
              } else if (scatter_type & ECollisionType::Scatter) {
                // update the throughput
                Float    pdf = tr_maj[heroChannel] * sigma_s[heroChannel];
                Spectrum throughput = tr_maj * sigma_s / pdf;
                beta *= throughput;
                r_u *= throughput;
                path.updateThroughput(throughput);

                if (rRec.depth++ >= m_maxDepth) {
                  terminated = true;

                  //* before terminate, connect to a subpath if rendering
                  if (!m_training) {
                    // TODO Actually a nee is also needed

                    L +=
                        beta * ReuseMultipleScattering(majRec.p) / (4.f * M_PI);
                  }

                  return false;
                }

                DirectSamplingRecord dRec{majRec.p, .0f};
                dRec.refN = Normal(.0f);

                Spectrum contrib = SampleVolumetricNEE(
                    scene, dRec, -ray.d, medium, rRec.sampler, r_u, heroChannel,
                    medium->getPhaseFunction(), nullptr);
                L += beta * contrib;
                path.addContribution(contrib);

                // sample a new direction
                const PhaseFunction        *phase = medium->getPhaseFunction();
                MediumSamplingRecord        mRec;
                PhaseFunctionSamplingRecord pRec{mRec, -ray.d};
                Float                       phasePdf;
                Float phaseWeight = phase->sample(pRec, phasePdf, rRec.sampler);

                Vector path_wi;
                Point  path_p;
                if (phaseWeight == .0f || phasePdf == .0f)
                  terminated = true;
                else {
                  scattered       = true;
                  specular_bounce = false;

                  beta *= phaseWeight;
                  path.updateThroughput(Spectrum(phaseWeight));

                  r_l      = r_u / phasePdf;
                  ray      = Ray{majRec.p, pRec.wo, .0f};
                  ray.mint = Epsilon;

                  prev_p = majRec.p;
                  prev_n = Normal(.0f);

                  path_wi = pRec.wo;
                  path_p  = majRec.p;
                }

                //* Add a new vertex into path
                path.addVertex(path_p, path_wi);

                return false;
              } else if (scatter_type & ECollisionType::Null) {
                //* A null scatter occurs, just keep tracking
                Spectrum sigma_maj = majRec.sigma_maj;

                Float    pdf = tr_maj[heroChannel] * sigma_n[heroChannel];
                Spectrum throughput = tr_maj * sigma_n / pdf;
                beta *= throughput;
                path.updateThroughput(throughput);

                r_u *= tr_maj * sigma_n / pdf;
                r_l *= tr_maj * sigma_maj / pdf;

                ray      = Ray{majRec.p, ray.d, .0f};
                ray.mint = Epsilon;

                return !beta.isZero() && !r_u.isZero();
              } else {
                Log(EError, "Shouldn't arrive here");
                exit(1);
              }
            });

        if (terminated) break;
        if (scattered) continue;

        Spectrum throughput = T_maj / T_maj[heroChannel];
        beta *= throughput;
        r_u *= throughput;
        r_l *= throughput;
        path.updateThroughput(throughput);
      }

      if (!its.isValid()) {
        Spectrum Le = scene->evalEnvironment(ray);
        if (rRec.depth == 1 || specular_bounce) {
          // First hit or pure specular bounce
          Spectrum contrib = beta * Le / r_u.average();
          L += beta * contrib;
          path.addContribution(contrib);
        } else {
          // Apply MIS
          Float p_l        = PDF_Nee(scene, prev_p, prev_n, ray);
          r_l              = r_l * p_l;
          Spectrum contrib = Le / (r_u + r_l).average();
          L += beta * contrib;
          path.addContribution(contrib);
        }
        break;
      }

      if (const Emitter *emitter = its.shape->getEmitter(); emitter) {
        Spectrum Le = emitter->eval(its, -ray.d);
        if (rRec.depth == 1 || specular_bounce) {
          // First hit or pure specular bounce
          Spectrum contrib = Le / r_u.average();
          L += beta * contrib;
          path.addContribution(contrib);
        } else {
          // Apply MIS
          Float p_l        = PDF_Nee(scene, prev_p, prev_n, ray, &its);
          r_l              = r_l * p_l;
          Spectrum contrib = Le / (r_u + r_l).average();
          L += beta * contrib;
          path.addContribution(contrib);
        }
      }

      const BSDF *bsdf = its.getBSDF(ray);
      if (bsdf->isNull()) {
        ray      = Ray{its.p + Epsilon * ray.d, ray.d, ray.time};
        ray.mint = Epsilon;
        if (its.isMediumTransition()) rRec.medium = its.getTargetMedium(ray.d);
        continue;
      }

      if (rRec.depth++ >= m_maxDepth) break;

      DirectSamplingRecord dRec(its);
      if (bsdf->getType() & BSDF::ESmooth) {
        Spectrum contrib =
            SampleVolumetricNEE(scene, dRec, -ray.d, rRec.medium, rRec.sampler,
                                r_u, heroChannel, nullptr, bsdf, &its);
        L += beta * contrib;
        path.addContribution(contrib);
      }

      // BSDF sampling
      BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
      Float              bsdfPdf;
      Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
      if (bsdfWeight.isZero()) break;

      const Vector wo        = its.toWorld(bRec.wo);
      Float        woDotGeoN = dot(its.geoFrame.n, wo);
      if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals) break;

      beta *= bsdfWeight;
      path.updateThroughput(bsdfWeight);

      ray = Ray(its.p + Epsilon * wo, wo, Epsilon, INFINITY, ray.time);
      specular_bounce  = (bRec.sampledType &
                         (BSDF::EDeltaReflection | BSDF::EDeltaTransmission));
      r_l              = r_u / bsdfPdf;
      prev_scatter_pdf = bsdfPdf;

      prev_p = its.p;
      prev_n = its.shFrame.n;

      if (its.isMediumTransition()) rRec.medium = its.getTargetMedium(ray.d);

      // TODO disable rr for simplicity
      //  if (rRec.depth >= m_rrDepth) {
      //    Float q = std::min(beta.max(), (Float)0.95f);
      //    if (rRec.nextSample1D() >= q) break;
      //    beta /= q;
      //  }
    }

    if (m_training) {
      // add the path vertex into a storage
      m_storage->addPath(path);
    }

    return L;
  }

  virtual bool preprocess(const Scene *scene, RenderQueue *queue,
                          const RenderJob *job, int sceneResID, int sensorResID,
                          int samplerResID) override;

protected:
  enum ECollisionType { Absorb = 0b0001, Scatter = 0b0010, Null = 0b0100 };

  ECollisionType SampleCollision(const Spectrum &sigma_a,
                                 const Spectrum &sigma_s,
                                 const Spectrum &sigma_n, Float u,
                                 uint32_t channel) const {
    // Only first channel here
    Float sum = sigma_a[channel] + sigma_s[channel] + sigma_n[channel];
    Float P_a = sigma_a[channel] / sum;
    Float P_s = sigma_s[channel] / sum;

    if (u < P_a) return ECollisionType::Absorb;
    if (u < P_a + P_s) return ECollisionType::Scatter;
    return ECollisionType::Null;
  }

  uint32_t SampleHeroChannel(Float u) const {
    return math::clamp((int)(u * 3), 0, 2);
  }

  Spectrum SampleVolumetricNEE(const Scene *scene, DirectSamplingRecord &dRec,
                               const Vector &wi, const Medium *medium,
                               Sampler *sampler, Spectrum r_p,
                               const uint32_t       heroChannel,
                               const PhaseFunction *phase = nullptr,
                               const BSDF          *bsdf  = nullptr,
                               const Intersection  *its   = nullptr) const {
    using RNG = RNG_CX;

    Spectrum Le{.0f};
    Spectrum scatterVal{.0f};
    Float    scatterPdf = .0f;
    Point2   sample     = sampler->next2D();

    Le = scene->sampleEmitterDirect(dRec, sample, false);

    if (Le.isZero() || dRec.pdf == .0f) return Spectrum(.0f);

    const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

    if (bsdf && (bsdf->getType() & BSDF::ESmooth)) {
      BSDFSamplingRecord bRec(*its, its->toLocal(dRec.d), ERadiance);
      scatterVal = bsdf->eval(bRec);
      scatterPdf = bsdf->pdf(bRec);
    } else if (phase) {
      MediumSamplingRecord        mRec; // TODO
      PhaseFunctionSamplingRecord pRec{mRec, wi, dRec.d};
      scatterVal = Spectrum{phase->eval(pRec)};
      scatterPdf = phase->pdf(pRec);
    }

    Spectrum        tr(1.f), r_l(1.f), r_u(1.f);
    RayDifferential shadowRay{dRec.ref, dRec.d, .0f};
    shadowRay.mint        = Epsilon;
    shadowRay.maxt        = dRec.dist * (1 - ShadowEpsilon);
    shadowRay.heroChannel = heroChannel;

    uint64_t      seed = ((uint64_t)rand() << 32) + (uint64_t)rand();
    RNG           rng{seed};
    const Medium *currentMedium = medium;
    Intersection  si;
    while (true) {

      scene->rayIntersect(shadowRay, si);

      if (si.isValid()) {
        if (const BSDF *bsdf = si.getBSDF(); !bsdf->isNull()) {
          // Hit opaque surface
          tr = Spectrum(.0f);
          break;
        }
      }

      if (currentMedium) {
        // Accumulate tr
        Float    tMax = si.isValid() ? si.t : shadowRay.maxt;
        Spectrum T_maj =
            SampleTrMajorant(currentMedium, shadowRay, tMax, rng,
                             ETrackingType::ENaiveDeltaTracking,
                             [&](MajorantSamplingRecord majRec) {
                               Spectrum sigma_n   = majRec.sigma_n;
                               Spectrum sigma_maj = majRec.sigma_maj;
                               Spectrum tr_maj    = majRec.tr_majorant;
                               Float    pdf =
                                   tr_maj[heroChannel] * sigma_maj[heroChannel];

                               tr *= tr_maj * sigma_n / pdf;
                               r_l *= tr_maj * sigma_maj / pdf;
                               r_u *= tr_maj * sigma_n / pdf;

                               if (tr.isZero()) return false;
                               return true;
                             });

        tr *= T_maj / T_maj[heroChannel];
        r_l *= T_maj / T_maj[heroChannel];
        r_u *= T_maj / T_maj[heroChannel];
      }

      if (tr.isZero()) break;

      shadowRay.o    = si.p;
      shadowRay.mint = Epsilon;
      shadowRay.maxt = shadowRay.maxt - si.t;

      if (!si.isValid()) break;

      if (si.isMediumTransition())
        currentMedium = si.getTargetMedium(shadowRay.d);
    }

    r_l *= r_p;
    r_u *= r_p * scatterPdf / dRec.pdf;

    if (!(emitter->isOnSurface() && dRec.measure == ESolidAngle))
      return scatterVal * tr * Le / r_l.average();
    else
      return scatterVal * tr * Le / (r_l + r_u).average();
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

  Spectrum ReuseMultipleScattering(Point p) const {
    //* nearest search
    const Float query_pt[3] = {p[0], p[1], p[2]};

    size_t                num_result = 16;
    std::vector<uint32_t> ret_idx(num_result);
    std::vector<Float>    out_dist_sqr(num_result);

    num_result = m_pathGrid->searchKNN(query_pt, num_result, ret_idx.data(),
                                       out_dist_sqr.data());

    Spectrum Ls(.0f);
    for (int i = 0; i < num_result; ++i) {
      auto [p_, v_, L_] = m_pathGrid->getData(ret_idx[0]);
      Ls += L_;
    }
    return Ls / num_result;
  }

private:
  //* Training variable
  int  m_training_spp;
  bool m_training;

  std::unique_ptr<pathgrid::PathGrid>    m_pathGrid;
  std::unique_ptr<pathgrid::PathStorage> m_storage;

  MTS_DECLARE_CLASS()
};

bool PathGrid::preprocess(const Scene *scene, RenderQueue *queue,
                          const RenderJob *job, int sceneResID, int sensorResID,
                          int samplerResID) {
  //* Generate a set of subpaths stored in kd-trees

  m_training = true;

  ref<Scheduler> sched = Scheduler::getInstance();
  ref<Sensor> sensor   = static_cast<Sensor *>(sched->getResource(sensorResID));
  ref<Scene>  train_scene =
      new Scene(static_cast<Scene *>(sched->getResource(sceneResID)));
  const int trainSceneResID = sched->registerResource(train_scene);

  // a loop ?
  {
    Properties training_sampler_props = scene->getSampler()->getProperties();

    training_sampler_props.removeProperty("sampleCount");
    training_sampler_props.setSize("sampleCount", m_training_spp);
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
    int training_sampler_resid = sched->registerMultiResource(samplers);
    for (size_t i = 0; i < sched->getCoreCount(); ++i)
      samplers[i]->decRef();

    SamplingIntegrator::render(train_scene, queue, job, trainSceneResID,
                               sensorResID, training_sampler_resid);

    sched->unregisterResource(training_sampler_resid);
    sensor->getFilm()->clear();

    m_pathGrid->addStorage(*m_storage);
    m_storage->clear();
  }

  m_pathGrid->init();

  Log(EInfo, "Preprocess end");
  m_training = false;

  return true;
}

MTS_IMPLEMENT_CLASS(PathGrid, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(PathGrid, "Volumetric PathGrid")

MTS_NAMESPACE_END