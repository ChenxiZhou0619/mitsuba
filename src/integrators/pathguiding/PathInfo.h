#pragma once
#include <mitsuba/mitsuba.h>
#include <openpgl/cpp/OpenPGL.h>
MTS_NAMESPACE_BEGIN

struct PathInfo {
  PathInfo(size_t max_size);

  void AddVertex(const Point3 &position);

  void AddLdContribution(const Spectrum &contribution);

  void UpdateBeta(const Spectrum &scatter_weight);

  void SetPdf(Float pdf);

  void SetDirection(const Vector3 &dir);

  void SetDistance(Float distance);

  size_t ToPGLSampleData(std::vector<PGLSampleData> &data) const;

  friend Spectrum FirstVertexLo(const PathInfo &p_info);

private:
  size_t size, max_size;

  std::vector<Point3>   positions;
  std::vector<Vector3>  directions;
  std::vector<Spectrum> contributions;
  std::vector<Float>    pdfs;
  std::vector<Float>    distances;
  std::vector<Spectrum> betas;

  // flags
};

inline PathInfo::PathInfo(size_t max_size) : size(0), max_size(max_size) {
  positions.reserve(max_size);
  directions.reserve(max_size);
  contributions.reserve(max_size);
  pdfs.reserve(max_size);
  distances.reserve(max_size);
  betas.reserve(max_size);
}

inline void PathInfo::AddVertex(const Point3 &position) {
  if (++size <= max_size) {
    positions.emplace_back(position);
    contributions.emplace_back(Spectrum(.0f));
    betas.emplace_back(Spectrum(1.f));
  }
}

inline void PathInfo::AddLdContribution(const Spectrum &contribution) {
  for (int idx = 0; idx < size - 1; ++idx) {
    contributions[idx] += contribution * betas[idx];
  }
}

inline void PathInfo::UpdateBeta(const Spectrum &beta) {
  for (int idx = 0; idx < size; ++idx) {
    betas[idx] = beta * betas[idx];

    //    betas[idx][0] = std::min(betas[idx][0], 10.f);
    //    betas[idx][1] = std::min(betas[idx][1], 10.f);
    //    betas[idx][2] = std::min(betas[idx][2], 10.f);
  }
}

inline void PathInfo::SetPdf(Float pdf) { pdfs.emplace_back(pdf); }

inline void PathInfo::SetDirection(const Vector3 &dir) { directions.emplace_back(dir); }

inline void PathInfo::SetDistance(Float distance) { distances.emplace_back(distance); }

inline size_t PathInfo::ToPGLSampleData(std::vector<PGLSampleData> &data) const {
  for (int idx = 0; idx < size; ++idx) {
    if (distances[idx] > 0 && !contributions[idx].isZero()) {
      PGLSampleData d;
      d.position  = {positions[idx][0], positions[idx][1], positions[idx][2]};
      d.direction = {directions[idx][0], directions[idx][1], directions[idx][2]};
      d.weight    = contributions[idx].average();
      d.pdf       = pdfs[idx];
      d.distance  = distances[idx];
      d.flags     = 0;
      data.emplace_back(d);
    }
  }
  return data.size();
}

inline Spectrum FirstVertexLo(const PathInfo &p_info) {
  return p_info.size > 0 ? p_info.contributions[0] : Spectrum(.0f);
}

MTS_NAMESPACE_END