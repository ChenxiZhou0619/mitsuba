/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

// clang-format off
#include <mitsuba/bidir/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/renderproc.h>
#include "mitsuba/core/plugin.h"

#include "gpt_proc.h"
#include "gpt_wr.h"
#include "../poisson_solver/Solver.hpp"


#include "gradient_pt.h"

#include "impl/pt.h"
#include "impl/pspt.h"
// clang-format on

MTS_NAMESPACE_BEGIN

/*!\plugin{gpt}{Gradient-domain path tracer}
* \order{5}
* \parameters{
*	   \parameter{reconstructL1}{\Boolean}{If set, the rendering method
reconstructs the final image using a reconstruction method
*           that efficiently kills many image artifacts. The reconstruction is
slightly biased, but the bias will go away by increasing sample count.
\default{\code{true}}
*     }
*	   \parameter{reconstructL2}{\Boolean}{If set, the rendering method
reconstructs the final image using a reconstruction method that is unbiased,
*			but sometimes introduces severe dipole artifacts.
\default{\code{false}}
*     }
*	   \parameter{shiftThreshold}{\Float}{Specifies the roughness threshold
for classifying materials as 'diffuse', in contrast to 'specular', *
for the purposes of constructing paths pairs for estimating pixel differences.
This value should usually be somewhere between 0.0005 and 0.01. *
If the result image has noise similar to standard path tracing, increasing or
decreasing this value may sometimes help. This implementation assumes that this
value is small.\default{\code{0.001}} *	   } *
\parameter{reconstructAlpha}{\Float}{ *			Higher value makes the
reconstruction trust the noisy color image more, giving less weight to the
usually lower-noise gradients.
*			The optimal value tends to be around 0.2, but for scenes
with much geometric detail at sub-pixel level a slightly higher value such as
0.3 or 0.4 may be tried.\default{\code{0.2}}
           }
* }
*
*
* This plugin implements a gradient-domain path tracer (short: G-PT) as
described in the paper "Gradient-Domain Path Tracing" by Kettunen et al.
* It samples difference images in addition to the standard color image, and
reconstructs the final image based on these.
* It supports classical materials like diffuse, specular and glossy materials,
and area and point lights, depth-of-field, and low discrepancy samplers.
* There is also experimental support for sub-surface scattering and motion blur.
Note that this is still an experimental implementation of Gradient-Domain Path
Tracing
* that has not been tested with all of Mitsuba's features. Notably there is no
support yet for any kind of participating media or directional lights.
* Environment maps are supported, though. Does not support the 'hide emitters'
option even though it is displayed.

*
*/

/// If defined, uses only the central sample for the throughput estimate.
/// Otherwise uses offset paths for estimating throughput too.
// #define CENTRAL_RADIANCE

/// If defined, applies reconstruction after rendering.
#define RECONSTRUCT

static StatsCounter avgPathLength("Gradient Path Tracer", "Average path length",
                                  EAverage);

// Output buffer names.
static const size_t BUFFER_FINAL =
    0; ///< Buffer index for the final image. Also used for preview.
static const size_t BUFFER_THROUGHPUT =
    1;                             ///< Buffer index for the noisy color image.
static const size_t BUFFER_DX = 2; ///< Buffer index for the X gradients.
static const size_t BUFFER_DY = 3; ///< Buffer index for the Y gradients.
static const size_t BUFFER_VERY_DIRECT =
    4; ///< Buffer index for very direct light.

GradientPathIntegrator::GradientPathIntegrator(const Properties &props)
    : MonteCarloIntegrator(props) {
  m_config.m_maxDepth       = props.getInteger("maxDepth", -1);
  m_config.m_minDepth       = props.getInteger("minDepth", -1);
  m_config.m_rrDepth        = props.getInteger("rrDepth", 5);
  m_config.m_strictNormals  = props.getBoolean("strictNormals", false);
  m_config.m_shiftThreshold = props.getFloat("shiftThreshold", Float(0.001));
  m_config.m_reconstructL1  = props.getBoolean("reconstructL1", true);
  m_config.m_reconstructL2  = props.getBoolean("reconstructL2", false);
  m_config.m_reconstructAlpha =
      (Float)props.getFloat("reconstructAlpha", Float(0.2));

  m_config.m_shiftType = props.getInteger("shiftType", 0);

  if (m_config.m_reconstructL1 && m_config.m_reconstructL2)
    Log(EError, "Disable 'reconstructL1' or 'reconstructL2': Cannot display "
                "two reconstructions at a time!");

  if (m_config.m_reconstructAlpha <= 0.0f)
    Log(EError, "'reconstructAlpha' must be set to a value greater than zero!");

  if (m_config.m_maxDepth <= 0 && m_config.m_maxDepth != -1)
    Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater "
                "than zero!");
}

GradientPathIntegrator::GradientPathIntegrator(Stream          *stream,
                                               InstanceManager *manager)
    : MonteCarloIntegrator(stream, manager) {
  m_config = GradientPathTracerConfig(stream);
}

void GradientPathIntegrator::renderBlock(
    const Scene *scene, const Sensor *sensor, Sampler *sampler,
    GPTWorkResult *block, const bool &stop,
    const std::vector<TPoint2<uint8_t>> &points) const {

  // TODO Choose different tracer here

  std::shared_ptr<IGradientPathTracer> tracer = nullptr;

  switch (m_config.m_shiftType) {
  case 0 /** Default half-vector + reconnect */:
    tracer = std::make_shared<GradientPathTracer>(scene, sensor, sampler, block,
                                                  &m_config);
    break;
  case 1 /** Primary space shift */:
    tracer = std::make_shared<PSGradientPathTracer>(scene, sensor, sampler,
                                                    block, &m_config);
    std::cout << "Primary space has not implemented!\n";
    exit(1);
    break;

  default:
    std::cout << "Fatal : Unknown shift type, terminate\n";
    exit(1);
  }

  bool needsApertureSample = sensor->needsApertureSample();
  bool needsTimeSample     = sensor->needsTimeSample();

  // Original code from SamplingIntegrator.
  Float diffScaleFactor = 1.0f / std::sqrt((Float)sampler->getSampleCount());

  // Get ready for sampling.
  RadianceQueryRecord rRec(scene, sampler);

  Point2          apertureSample(0.5f);
  Float           timeSample = 0.5f;
  RayDifferential sensorRay;

  block->clear();

  // Sample at the given positions.
  Spectrum gradients[4];
  Spectrum shiftedThroughputs[4];

  for (size_t i = 0; i < points.size(); ++i) {
    if (stop) {
      break;
    }

    Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
    sampler->generate(offset);

    for (size_t j = 0; j < sampler->getSampleCount(); ++j) {
      if (stop) {
        break;
      }

      // Get the initial ray to sample.
      rRec.newQuery(RadianceQueryRecord::ESensorRay, sensor->getMedium());

      Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

      if (needsApertureSample) {
        apertureSample = rRec.nextSample2D();
      }
      if (needsTimeSample) {
        timeSample = rRec.nextSample1D();
      }

      // Do the actual sampling.
      Spectrum centralVeryDirect = Spectrum(0.0f);
      Spectrum centralThroughput = Spectrum(0.0f);

      tracer->evaluatePoint(rRec, samplePos, apertureSample, timeSample,
                            diffScaleFactor, centralVeryDirect,
                            centralThroughput, gradients, shiftedThroughputs);

      // Accumulate results.
      const Point2 right_pixel  = samplePos + Vector2(1.0f, 0.0f);
      const Point2 bottom_pixel = samplePos + Vector2(0.0f, 1.0f);
      const Point2 left_pixel   = samplePos - Vector2(1.0f, 0.0f);
      const Point2 top_pixel    = samplePos - Vector2(0.0f, 1.0f);
      const Point2 center_pixel = samplePos;

      static const int RIGHT  = 0;
      static const int BOTTOM = 1;
      static const int LEFT   = 2;
      static const int TOP    = 3;

      // Note: Sampling differences and throughputs to multiple directions is
      // essentially
      //       multiple importance sampling (MIS) between the pixels.
      //
      //       For a sample from a strategy participating in the MIS to be
      //       unbiased, we need to divide its result by the selection
      //       probability of that strategy.
      //
      //       As the selection probability is 0.5 for both directions (no
      //       adaptive sampling), we need to multiply the results by two.

      // Note: The central pixel is estimated as
      //               1/4 * (throughput estimate sampled from MIS(center, top)
      //                      + throughput estimate sampled from MIS(center,
      //                      right)
      //                      + throughput estimate sampled from MIS(center,
      //                      bottom)
      //                      + throughput estimate sampled from MIS(center,
      //                      left)).
      //
      //       Variable centralThroughput is the sum of four throughput
      //       estimates sampled from each of these distributions, from the
      //       central pixel, so it's actually four samples, and thus its weight
      //       is 4.
      //
      //       The other samples from the MIS'd distributions will be sampled
      //       from the neighboring pixels, and their weight is 1.
      //
      //       If this feels too complicated, it should be OK to output a
      //       standard throughput sample from the path tracer.

      // Add the throughput image as a preview. Note: Preview and final buffers
      // are shared.
      {
#ifdef CENTRAL_RADIANCE
        block->put(samplePos, centralVeryDirect + centralThroughput, 1.0f, 1.0f,
                   BUFFER_FINAL); // Standard throughput estimate with direct.
#else
        block->put(samplePos, (8 * centralVeryDirect) + (2 * centralThroughput),
                   4.0f, 4.0f,
                   0); // Adds very direct on top of the throughput image.

        block->put(left_pixel, (2 * shiftedThroughputs[LEFT]), 1.0f, 1.0f,
                   BUFFER_FINAL); // Negative x throughput.
        block->put(right_pixel, (2 * shiftedThroughputs[RIGHT]), 1.0f, 1.0f,
                   BUFFER_FINAL); // Positive x throughput.
        block->put(top_pixel, (2 * shiftedThroughputs[TOP]), 1.0f, 1.0f,
                   BUFFER_FINAL); // Negative y throughput.
        block->put(bottom_pixel, (2 * shiftedThroughputs[BOTTOM]), 1.0f, 1.0f,
                   BUFFER_FINAL); // Positive y throughput.
#endif
      }

      // Actual throughputs, with MIS between central and neighbor pixels for
      // all neighbors. This can be replaced with a standard throughput sample
      // without much loss of quality in most cases.
      {
#ifdef CENTRAL_RADIANCE
        block->put(samplePos, centralThroughput, 1.0f, 1.0f,
                   BUFFER_THROUGHPUT); // Standard throughput estimate.
#else
        block->put(samplePos, (2 * centralThroughput), 4.0f, 4.0f,
                   BUFFER_THROUGHPUT); // Central throughput.

        block->put(left_pixel, (2 * shiftedThroughputs[LEFT]), 1.0f, 1.0f,
                   BUFFER_THROUGHPUT); // Negative x throughput.
        block->put(right_pixel, (2 * shiftedThroughputs[RIGHT]), 1.0f, 1.0f,
                   BUFFER_THROUGHPUT); // Positive x throughput.
        block->put(top_pixel, (2 * shiftedThroughputs[TOP]), 1.0f, 1.0f,
                   BUFFER_THROUGHPUT); // Negative y throughput.
        block->put(bottom_pixel, (2 * shiftedThroughputs[BOTTOM]), 1.0f, 1.0f,
                   BUFFER_THROUGHPUT); // Positive y throughput.
#endif
      }

      // Gradients.
      {
        block->put(left_pixel, -(2 * gradients[LEFT]), 1.0f, 1.0f,
                   BUFFER_DX); // Negative x gradient.
        block->put(center_pixel, (2 * gradients[RIGHT]), 1.0f, 1.0f,
                   BUFFER_DX); // Positive x gradient.
        block->put(top_pixel, -(2 * gradients[TOP]), 1.0f, 1.0f,
                   BUFFER_DY); // Negative y gradient.
        block->put(center_pixel, (2 * gradients[BOTTOM]), 1.0f, 1.0f,
                   BUFFER_DY); // Positive y gradient.
      }

      // Very direct.
      block->put(center_pixel, centralVeryDirect, 1.0f, 1.0f,
                 BUFFER_VERY_DIRECT);
    }
  }
}

/// Custom render function that samples a number of paths for evaluating
/// differences between pixels.
bool GradientPathIntegrator::render(Scene *scene, RenderQueue *queue,
                                    const RenderJob *job, int sceneResID,
                                    int sensorResID, int samplerResID) {
  if (m_hideEmitters) {
    /* Not supported! */
    Log(EError, "Option 'hideEmitters' not implemented for Gradient-Domain "
                "Path Tracing!");
  }

  /* Get config from the parent class. */
  m_config.m_maxDepth      = m_maxDepth;
  m_config.m_minDepth      = 1; // m_minDepth;
  m_config.m_rrDepth       = m_rrDepth;
  m_config.m_strictNormals = m_strictNormals;

  /* Code duplicated from SamplingIntegrator::Render. */
  ref<Scheduler> sched = Scheduler::getInstance();
  ref<Sensor> sensor   = static_cast<Sensor *>(sched->getResource(sensorResID));

  /* Set up MultiFilm. */
  ref<Film> film = sensor->getFilm();

  std::vector<std::string> outNames = {"-final", "-throughput", "-dx", "-dy",
                                       "-direct"};
  if (!film->setBuffers(outNames)) {
    Log(EError, "Cannot render image! G-PT has been called without MultiFilm.");
    return false;
  }

  size_t         nCores = sched->getCoreCount();
  const Sampler *sampler =
      static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
  size_t sampleCount = sampler->getSampleCount();

  Log(EInfo,
      "Starting render job (GPT::render) (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
      " %s, " SSE_STR ") ..",
      film->getCropSize().x, film->getCropSize().y, sampleCount,
      sampleCount == 1 ? "sample" : "samples", nCores,
      nCores == 1 ? "core" : "cores");

  /* This is a sampling-based integrator - parallelize. */
  ref<BlockedRenderProcess> proc =
      new GPTRenderProcess(job, queue, scene->getBlockSize(), m_config);

  int integratorResID = sched->registerResource(this);
  proc->bindResource("integrator", integratorResID);
  proc->bindResource("scene", sceneResID);
  proc->bindResource("sensor", sensorResID);
  proc->bindResource("sampler", samplerResID);

  scene->bindUsedResources(proc);
  bindUsedResources(proc);
  sched->schedule(proc);

  m_process = proc;
  sched->wait(proc);

  sched->unregisterResource(integratorResID);
  m_process = NULL;

#ifdef RECONSTRUCT
  /* Reconstruct. */
  if (m_config.m_reconstructL1 || m_config.m_reconstructL2) {
    /* Allocate bitmaps for the solvers. */
    ref<Bitmap> throughputBitmap(
        new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));
    ref<Bitmap> directBitmap(
        new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));
    ref<Bitmap> dxBitmap(
        new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));
    ref<Bitmap> dyBitmap(
        new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));
    ref<Bitmap> reconstructionBitmap(
        new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));

    /* Develop primal and gradient data into bitmaps. */
    film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0),
                       throughputBitmap, BUFFER_THROUGHPUT);
    film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0),
                       dxBitmap, BUFFER_DX);
    film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0),
                       dyBitmap, BUFFER_DY);
    film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0),
                       directBitmap, BUFFER_VERY_DIRECT);

    /* Transform the data for the solver. */
    size_t subPixelCount = 3 * film->getCropSize().x * film->getCropSize().y;
    std::vector<float> throughputVector(subPixelCount);
    std::vector<float> dxVector(subPixelCount);
    std::vector<float> dyVector(subPixelCount);
    std::vector<float> directVector(subPixelCount);
    std::vector<float> reconstructionVector(subPixelCount);

    std::transform(throughputBitmap->getFloatData(),
                   throughputBitmap->getFloatData() + subPixelCount,
                   throughputVector.begin(), [](Float x) { return (float)x; });
    std::transform(dxBitmap->getFloatData(),
                   dxBitmap->getFloatData() + subPixelCount, dxVector.begin(),
                   [](Float x) { return (float)x; });
    std::transform(dyBitmap->getFloatData(),
                   dyBitmap->getFloatData() + subPixelCount, dyVector.begin(),
                   [](Float x) { return (float)x; });
    std::transform(directBitmap->getFloatData(),
                   directBitmap->getFloatData() + subPixelCount,
                   directVector.begin(), [](Float x) { return (float)x; });

    /* Reconstruct. */
    poisson::Solver::Params params;

    if (m_config.m_reconstructL1) {
      params.setConfigPreset("L1D");
    } else if (m_config.m_reconstructL2) {
      params.setConfigPreset("L2D");
    }

    params.alpha = (float)m_config.m_reconstructAlpha;
    params.setLogFunction(
        poisson::Solver::Params::LogFunction([](const std::string &message) {
          SLog(EInfo, "%s", message.c_str());
        }));

    poisson::Solver solver(params);
    solver.importImagesMTS(dxVector.data(), dyVector.data(),
                           throughputVector.data(), directVector.data(),
                           film->getCropSize().x, film->getCropSize().y);

    solver.setupBackend();
    solver.solveIndirect();

    solver.exportImagesMTS(reconstructionVector.data());

    /* Give the solution back to Mitsuba. */
    int w = reconstructionBitmap->getSize().x;
    int h = reconstructionBitmap->getSize().y;

    for (int y = 0, p = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x, p += 3) {
        Float color[3] = {(Float)reconstructionVector[p],
                          (Float)reconstructionVector[p + 1],
                          (Float)reconstructionVector[p + 2]};
        reconstructionBitmap->setPixel(Point2i(x, y), Spectrum(color));
      }
    }

    film->setBitmapMulti(reconstructionBitmap, 1, BUFFER_FINAL);
  }
#endif

  return proc->getReturnStatus() == ParallelProcess::ESuccess;
}

static Float miWeight(Float pdfA, Float pdfB) {
  pdfA *= pdfA;
  pdfB *= pdfB;
  return pdfA / (pdfA + pdfB);
}

Spectrum GradientPathIntegrator::Li(const RayDifferential &r,
                                    RadianceQueryRecord   &rRec) const {
  // Duplicate of MIPathTracer::Li to support sub-surface scattering
  // initialization.

  /* Some aliases and local variables */
  const Scene    *scene = rRec.scene;
  Intersection   &its   = rRec.its;
  RayDifferential ray(r);
  Spectrum        Li(0.0f);
  bool            scattered = false;

  /* Perform the first ray intersection (or ignore if the
          intersection has already been provided). */
  rRec.rayIntersect(ray);
  ray.mint = Epsilon;

  Spectrum throughput(1.0f);
  Float    eta = 1.0f;

  while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
    if (!its.isValid()) {
      /* If no intersection could be found, potentially return
              radiance from a environment luminaire if it exists */
      if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) &&
          (!m_hideEmitters || scattered))
        Li += throughput * scene->evalEnvironment(ray);
      break;
    }

    const BSDF *bsdf = its.getBSDF(ray);

    /* Possibly include emitted radiance if requested */
    if (its.isEmitter() &&
        (rRec.type & RadianceQueryRecord::EEmittedRadiance) &&
        (!m_hideEmitters || scattered))
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
    Float              bsdfPdf;
    BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
    Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
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

    /* Keep track of the throughput and relative
            refractive index along the path */
    throughput *= bsdfWeight;
    eta *= bRec.eta;

    /* If a luminaire was hit, estimate the local illumination and
            weight using the power heuristic */
    if (hitEmitter &&
        (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
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
      if (rRec.nextSample1D() >= q) break;
      throughput /= q;
    }
  }

  return Li;
}

void GradientPathIntegrator::serialize(Stream          *stream,
                                       InstanceManager *manager) const {
  MonteCarloIntegrator::serialize(stream, manager);
  m_config.serialize(stream);
}

std::string GradientPathIntegrator::toString() const {
  std::ostringstream oss;
  oss << "GradientPathTracer[" << endl
      << "  maxDepth = " << m_maxDepth << "," << endl
      << "  rrDepth = " << m_rrDepth << "," << endl
      << "  shiftThreshold = " << m_config.m_shiftThreshold << endl
      << "  reconstructL1 = " << m_config.m_reconstructL1 << endl
      << "  reconstuctL2 = " << m_config.m_reconstructL2 << endl
      << "  reconstructAlpha = " << m_config.m_reconstructAlpha << endl
      << "]";
  return oss.str();
}

MTS_IMPLEMENT_CLASS_S(GradientPathIntegrator, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(GradientPathIntegrator, "Gradient Path Integrator");
MTS_NAMESPACE_END
