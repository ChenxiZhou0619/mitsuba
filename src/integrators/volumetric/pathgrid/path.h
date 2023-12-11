// clang-format off
#pragma once

#include <mitsuba/mitsuba.h>
#include <mutex>
// clang-format on

MTS_NAMESPACE_BEGIN

namespace pathgrid {

struct PathInfo {
  std::vector<Point>    ps;
  std::vector<Vector>   wis;
  std::vector<Spectrum> contribs;
  std::vector<Spectrum> betas;

  int max_size;
  int size;

  PathInfo(int max_size) : max_size(max_size), size(0) {
    ps.reserve(max_size);
    wis.reserve(max_size);
    contribs.reserve(max_size);
    betas.reserve(max_size);
  }

  void addVertex(Point p, Vector wi) {
    size++;
    ps.emplace_back(p);
    wis.emplace_back(wi);
    contribs.emplace_back(Spectrum(.0f));
    betas.emplace_back(Spectrum(1.f));
  }

  void addContribution(Spectrum contrib) {
    for (int i = 0; i < size; ++i) {
      contribs[i] += betas[i] * contrib;
    }
  }

  void updateThroughput(Spectrum beta) {
    for (int i = 0; i < size; ++i) {
      betas[i] *= beta;
    }
  }
};

struct PathStorage {
public:
  PathStorage(int max_size) : max_size(max_size), size(0) {
    ps.reserve(max_size);
    wis.reserve(max_size);
    contribs.reserve(max_size);
  }

  bool addPath(const PathInfo &path) {
    std::lock_guard<std::mutex> lock(mtx);

    for (int i = 1; i < path.size; ++i) {
      if (path.contribs[i].getLuminance() < .1f) continue;

      if (size++ >= max_size) return false;
      ps.emplace_back(path.ps[i]);
      wis.emplace_back(path.wis[i]);
      contribs.emplace_back(path.contribs[i]);
    }
    return true;
  }

  void clear() {
    size = 0;
    ps.clear();
    wis.clear();
    contribs.clear();
  }

public:
  int        max_size, size;
  std::mutex mtx;

  std::vector<Point>    ps;
  std::vector<Vector>   wis;
  std::vector<Spectrum> contribs;
};

} // namespace pathgrid

MTS_NAMESPACE_END