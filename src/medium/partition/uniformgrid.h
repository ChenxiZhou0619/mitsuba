#pragma once

#include "partition.h"
#include <tbb/tbb.h>
MTS_NAMESPACE_BEGIN

class UniformGridTracker : public Tracker {
public:
  UniformGridTracker() = default;

  UniformGridTracker(const Vector3i &resolution, const Vector &voxelSize,
                     const Point &startVoxel, const Vector &direction, Float t,
                     Float tmax);

  virtual std::optional<SegmentData> nextSeg() override;

  virtual void marchForward(Float step) override;

  virtual void marchNext() override;

private:
  bool  terminate;
  Float tmax;
  Float nextCrossingT[3];
  Float deltaT[3];
  int   step[3];
  int   voxelLimit[3];
  int   currentIdx[3];
  int   stepAxis;

  Vector3i resolution;
};

class UniformGrid : public PartitionGrid {
public:
  UniformGrid() = delete;

  UniformGrid(const AABB &rootBound, Vector3i resolution)
      : PartitionGrid(rootBound), m_resolution(resolution) {

    Float x = (rootBound.max[0] - rootBound.min[0]) / resolution[0];
    Float y = (rootBound.max[1] - rootBound.min[1]) / resolution[1];
    Float z = (rootBound.max[2] - rootBound.min[2]) / resolution[2];

    m_voxelSize = Vector{x, y, z};
    m_majData =
        std::vector<Float>(resolution[0] * resolution[1] * resolution[2]);
  }

  virtual void
  configure(const std::function<Point3(Point)>  &w2i,
            const std::function<Float(Point3i)> &accessor) override;

  virtual void accessData(MajorantInfo *info,
                          uint32_t      dataAccessIdx) const override;

  virtual std::unique_ptr<Tracker> getTracker(const Ray &ray,
                                              Float      tmax) const override;

protected:
  AABB getVoxelBound(int voxelIndex) const {
    int z = voxelIndex % m_resolution[2];
    int y = (voxelIndex / m_resolution[2]) % m_resolution[1];
    int x = voxelIndex / (m_resolution[1] * m_resolution[2]);

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
                               m_rootBound.min, m_rootBound.max);
    Point voxelBoundMax = Lerp({xRatioMax, yRatioMax, zRatioMax},
                               m_rootBound.min, m_rootBound.max);

    return AABB{voxelBoundMin, voxelBoundMax};
  }

private:
  Vector3i           m_resolution;
  Vector3            m_voxelSize;
  std::vector<Float> m_majData;
};

inline UniformGridTracker::UniformGridTracker(const Vector3i &resolution,
                                              const Vector   &voxelSize,
                                              const Point    &startVoxel,
                                              const Vector &direction, Float t,
                                              Float tmax) {
  this->t          = t;
  this->tmax       = tmax;
  this->resolution = resolution;
  terminate        = false;
  stepAxis         = -1;

  // Initialization
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

inline std::optional<SegmentData> UniformGridTracker::nextSeg() {
  auto voxelToOffset = [&](int x, int y, int z) {
    int offset = x * (resolution[1] * resolution[2]) + y * resolution[2] + z;
    return offset;
  };

  if (terminate) return std::nullopt;

  if (nextCrossingT[0] < nextCrossingT[1] &&
      nextCrossingT[0] < nextCrossingT[2])
    stepAxis = 0;
  else if (nextCrossingT[1] < nextCrossingT[2])
    stepAxis = 1;
  else
    stepAxis = 2;

  SegmentData seg;
  seg.dataAccessIdx =
      voxelToOffset(currentIdx[0], currentIdx[1], currentIdx[2]);
  seg.segDistance = nextCrossingT[stepAxis] - t;

  if (nextCrossingT[stepAxis] > tmax) {
    seg.segDistance = tmax - t;
  }

  return std::make_optional(seg);
}

inline void UniformGridTracker::marchForward(Float step) {
  t += step;
  if (t >= tmax) terminate = true;
}

inline void UniformGridTracker::marchNext() {
  t = nextCrossingT[stepAxis];
  currentIdx[stepAxis] += step[stepAxis];
  nextCrossingT[stepAxis] += deltaT[stepAxis];

  if (currentIdx[stepAxis] == voxelLimit[stepAxis] || t > tmax)
    terminate = true;

  stepAxis = -1;
}

inline void
UniformGrid::configure(const std::function<Point3(Point)>  &w2i,
                       const std::function<Float(Point3i)> &accessor) {

  int nVoxels = m_resolution[0] * m_resolution[1] * m_resolution[2];

  tbb::parallel_for(
      tbb::blocked_range<int>(0, nVoxels), [&](tbb::blocked_range<int> range) {
        for (int r = range.begin(); r < range.end(); ++r) {

          auto [pmin, pmax] = getVoxelBound(r);
          Float x0 = INFINITY, y0 = INFINITY, z0 = INFINITY;
          Float x1 = -INFINITY, y1 = -INFINITY, z1 = -INFINITY;

          for (int i = 0; i < 8; ++i) {
            Point pvtx;
            pvtx[0] = (i & 0b0001) ? pmax[0] : pmin[0];
            pvtx[1] = (i & 0b0010) ? pmax[1] : pmin[1];
            pvtx[2] = (i & 0b0100) ? pmax[2] : pmin[2];

            pvtx = w2i(pvtx);

            x0 = std::min(x0, pvtx[0]);
            y0 = std::min(y0, pvtx[1]);
            z0 = std::min(z0, pvtx[2]);

            x1 = std::max(x1, pvtx[0]);
            y1 = std::max(y1, pvtx[1]);
            z1 = std::max(z1, pvtx[2]);
          }

          int ix0 = (int)(x0 - 1);
          int iy0 = (int)(y0 - 1);
          int iz0 = (int)(z0 - 1);

          int ix1 = (int)(x1 + 1);
          int iy1 = (int)(y1 + 1);
          int iz1 = (int)(z1 + 1);

          //          int ix0 = x0;
          //          int iy0 = y0;
          //          int iz0 = z0;
          //
          //          int ix1 = x1 + 1;
          //          int iy1 = y1 + 1;
          //          int iz1 = z1 + 1;

          Float majorant = .0f;
          for (int xx = ix0; xx <= ix1; ++xx)
            for (int yy = iy0; yy <= iy1; ++yy)
              for (int zz = iz0; zz <= iz1; ++zz) {
                majorant = std::max(majorant, accessor({xx, yy, zz}));
              }

          m_majData[r] = majorant;
        }
      });
}

inline void UniformGrid::accessData(MajorantInfo *info,
                                    uint32_t      dataAccessIdx) const {
  info->majorantDensity = m_majData[dataAccessIdx];
}

inline std::unique_ptr<Tracker> UniformGrid::getTracker(const Ray &ray,
                                                        Float      tmax) const {
  Float tnear, tfar;
  bool  overlap = m_rootBound.rayIntersect(ray, tnear, tfar);
  if (!overlap || tfar < 0) return nullptr;

  Float tmin = tnear > 0 ? tnear : .0f;
  tmax       = std::min(tmax, tfar);

  auto worldToIndex = [&](Point world) {
    Float xOffset = world[0] - m_rootBound.min[0];
    Float yOffset = world[1] - m_rootBound.min[1];
    Float zOffset = world[2] - m_rootBound.min[2];

    return Point3{xOffset / m_voxelSize[0], yOffset / m_voxelSize[1],
                  zOffset / m_voxelSize[2]};
  };

  Point startVoxel = worldToIndex(ray(tmin));

  // Initialize the UniformGridTracker
  return std::make_unique<UniformGridTracker>(m_resolution, m_voxelSize,
                                              startVoxel, ray.d, tmin, tmax);
}

MTS_NAMESPACE_END