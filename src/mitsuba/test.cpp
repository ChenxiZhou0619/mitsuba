#include <matplotlibcpp.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/mitsuba.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/render/scenehandler.h>

#include "../medium/partition/partition.h"
#include <mitsuba/render/medium.h>

namespace plt = matplotlibcpp;
using namespace mitsuba;

int main(int argc, char **argv) {

  //* System initialization
  Class::staticInitialization();
  Object::staticInitialization();
  PluginManager::staticInitialization();
  Statistics::staticInitialization();
  Thread::staticInitialization();
  Logger::staticInitialization();
  FileStream::staticInitialization();
  Spectrum::staticInitialization();
  Bitmap::staticInitialization();
  Scheduler::staticInitialization();
  SHVector::staticInitialization();
  SceneHandler::staticInitialization();

  PluginManager *manager = PluginManager::getInstance();

  //* Hard code medium configuration
  Properties p;
  p.setPluginName("nanovdbmedium");
  p.setTransform("toWorld", Transform());
  p.setString("filename", "/home/chenxizhou/Desktop/Programming/mitsuba/scenes/"
                          "disney-cloud/models/cloud.nvdb");
  p.setSpectrum("sigmaS", Spectrum(1.f));
  p.setSpectrum("sigmaA", Spectrum(0.f));

  constexpr int resolution = 128;
  p.setInteger("majResolutionX", resolution);
  p.setInteger("majResolutionY", resolution);
  p.setInteger("majResolutionZ", resolution);

  p.setString("partition", "octreeGrid");

  auto medium = (Medium *)manager->createObject(p);

  medium->configure();

  Ray    testRay;
  Point  origin    = Point(709.657, 382.659, 29.2989);
  Vector direction = normalize(
      Vector(-0.8930000000000291, -0.4479999999999791, -0.03399999999999892));

  testRay.setOrigin(origin);
  testRay.setDirection(direction);

  auto tracker = medium->GetTracker(testRay, INFINITY);

  std::vector<Float> t;
  std::vector<Float> maj;

  int prev_segvoxel = -1;

  // Record the majorant information
  while (auto segOpt = tracker->nextSeg()) {
    auto [segDistance, segVoxel] = *segOpt;

    Float segMin = tracker->t;
    Float val    = medium->SampleMajDensity(segVoxel);

    t.emplace_back(segMin);
    maj.emplace_back(val);

    t.emplace_back(segMin + segDistance);
    maj.emplace_back(val);

    tracker->marchNext();
  }

  Float tmin = t[0];
  Float tmax = t.back();

  std::vector<Float> t_prime;
  std::vector<Float> density;
  constexpr float    delta = .5f;
  while (tmin < tmax) {
    t_prime.emplace_back(tmin);

    Point p = testRay(tmin);
    density.emplace_back(medium->SampleDensity(p));

    tmin += delta;
  }

  plt::figure_size(12000, 720);
  plt::named_plot("majorant", t, maj);
  plt::named_plot("density", t_prime, density);

  plt::title("Partition along ray");
  plt::legend();
  plt::save("./maj.png");
}
