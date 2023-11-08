#include "majorantgrid.h"
#include <mitsuba/core/aabb.h>
MTS_NAMESPACE_BEGIN

MajorantGrid::MajorantGrid(const AABB     &mediumWorldBound,
                           const Vector3i &resolution) {
  m_resolution       = resolution;
  m_mediumWorldBound = mediumWorldBound;

  Float xSize =
      (mediumWorldBound.max[0] - mediumWorldBound.min[0]) / resolution[0];
  Float ySize =
      (mediumWorldBound.max[1] - mediumWorldBound.min[1]) / resolution[1];
  Float zSize =
      (mediumWorldBound.max[2] - mediumWorldBound.min[2]) / resolution[2];

  m_voxelSizew = {xSize, ySize, zSize};

  m_majData =
      std::vector<Float>(m_resolution[0] * m_resolution[1] * m_resolution[2]);
}

Float &MajorantGrid::at(int x, int y, int z) {
  // x * (resolution.y * resolution.z) + y * resolution.z + z
  int offset =
      x * (m_resolution[1] * m_resolution[2]) + y * m_resolution[2] + z;
  return m_majData[offset];
}

Float MajorantGrid::at(int x, int y, int z) const {
  // x * (resolution.y * resolution.z) + y * resolution.z + z
  int offset =
      x * (m_resolution[1] * m_resolution[2]) + y * m_resolution[2] + z;
  return m_majData[offset];
}

Vector3i MajorantGrid::offsetToXYZ(int offset) const {
  int z = offset % m_resolution[2];
  int y = (offset / m_resolution[2]) % m_resolution[1];
  int x = offset / (m_resolution[1] * m_resolution[2]);
  return {x, y, z};
}

AABB MajorantGrid::voxelBound(int x, int y, int z) const {
  auto Lerp = [](Point ratio, Point p1, Point p2) -> Point {
    Float x = (1 - ratio[0]) * p1[0] + ratio[0] * p2[0];
    Float y = (1 - ratio[1]) * p1[1] + ratio[1] * p2[1];
    Float z = (1 - ratio[2]) * p1[2] + ratio[2] * p2[2];

    return {x, y, z};
  };

  Float xRatioMin = (Float)x / m_resolution[0];
  Float yRatioMin = (Float)y / m_resolution[0];
  Float zRatioMin = (Float)z / m_resolution[0];

  Float xRatioMax = (Float)(x + 1) / m_resolution[0];
  Float yRatioMax = (Float)(y + 1) / m_resolution[0];
  Float zRatioMax = (Float)(z + 1) / m_resolution[0];

  Point voxelBoundMin = Lerp({xRatioMin, yRatioMin, zRatioMin},
                             m_mediumWorldBound.min, m_mediumWorldBound.max);
  Point voxelBoundMax = Lerp({xRatioMax, yRatioMax, zRatioMax},
                             m_mediumWorldBound.min, m_mediumWorldBound.max);

  return AABB{voxelBoundMin, voxelBoundMax};
}

Point MajorantGrid::worldToIndex(Point world) const {
  Float xOffset = world[0] - m_mediumWorldBound.min[0];
  Float yOffset = world[1] - m_mediumWorldBound.min[1];
  Float zOffset = world[2] - m_mediumWorldBound.min[2];

  return {xOffset / m_voxelSizew[0], yOffset / m_voxelSizew[1],
          zOffset / m_voxelSizew[2]};
}

DDATracker::DDATracker(const Vector3i &resolution, const Vector &voxelSize,
                       const Point &startVoxel, const Vector &direction,
                       Float t, Float tmax)
    : t(t), terminate(false), tmax(tmax), stepAxis(-1) {
  for (int axis = 0; axis < 3; ++axis) {
    currentIdx[axis] =
        math::clamp((int)startVoxel[axis], 0, resolution[axis] - 1);
    deltaT[axis] = voxelSize[axis] / std::abs(direction[axis]);

    if (direction[axis] >= .0f) {
      nextCrossingT[axis] = t + (currentIdx[axis] + 1.f - startVoxel[axis]) *
                                    voxelSize[axis] / direction[axis];
      step[axis]       = 1;
      voxelLimit[axis] = resolution[axis];
    } else {
      nextCrossingT[axis] = t + (currentIdx[axis] - startVoxel[axis]) *
                                    voxelSize[axis] / direction[axis];
      step[axis]       = -1;
      voxelLimit[axis] = -1;
    }
  }
}

std::optional<DDASegement> DDATracker::nextSeg() {
  if (terminate) return std::nullopt;

  if (nextCrossingT[0] < nextCrossingT[1] &&
      nextCrossingT[0] < nextCrossingT[2])
    stepAxis = 0;
  else if (nextCrossingT[1] < nextCrossingT[2])
    stepAxis = 1;
  else
    stepAxis = 2;

  DDASegement seg;
  seg.segVoxel    = {currentIdx[0], currentIdx[1], currentIdx[2]};
  seg.segDistance = nextCrossingT[stepAxis] - t;

  if (nextCrossingT[stepAxis] > tmax) {
    seg.segDistance = tmax - t;
  }

  return std::make_optional(seg);
}

void DDATracker::marchForward(Float step) {
  // Just march t, no update for voxel
  t += step;
  if (t >= tmax) terminate = true;
}

void DDATracker::marchNext() {

  t = nextCrossingT[stepAxis];
  currentIdx[stepAxis] += step[stepAxis];
  nextCrossingT[stepAxis] += deltaT[stepAxis];

  if (currentIdx[stepAxis] == voxelLimit[stepAxis] || t > tmax)
    terminate = true;

  stepAxis = -1;
}

MTS_NAMESPACE_END