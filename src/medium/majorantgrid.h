#pragma once
#include <mitsuba/core/aabb.h>
#include <mitsuba/mitsuba.h>
#include <optional>
MTS_NAMESPACE_BEGIN
//* A worldspace axis-aligned unified voxel grid
struct MajorantGrid {
  MajorantGrid() = delete;

  MajorantGrid(const AABB &mediumWorldBound, const Vector3i &resolution);

  Float &at(int x, int y, int z);

  Float at(int x, int y, int z) const;

  size_t size() const { return m_majData.size(); }

  Vector3i offsetToXYZ(int offset) const;

  AABB voxelBound(int x, int y, int z) const;

  Point worldToIndex(Point world) const;

  AABB               m_mediumWorldBound;
  Vector3i           m_resolution; // resolution of the majGrid
  Vector             m_voxelSizew; // size of a single voxel
  std::vector<Float> m_majData;    // store majorants
};

struct DDASegement {
  Float    segDistance;
  Vector3i segVoxel;
};

struct DDATracker {
  DDATracker() : terminate(true) {}

  DDATracker(const Vector3i &resolution, const Vector &voxelSize,
             const Point &startVoxel, const Vector &direction, Float t,
             Float tmax);

  std::optional<DDASegement> nextSeg();

  void marchForward(Float step);

  void marchNext();

  Float t;

private:
  bool  terminate;
  Float tmax;
  Float nextCrossingT[3];
  Float deltaT[3];
  int   step[3];
  int   voxelLimit[3];
  int   currentIdx[3];
  int   stepAxis;
};

MTS_NAMESPACE_END