/**
 *  Density spatitial partition structure, reduce the null collision in each
 * partitioned space
 */
#pragma once
#include <functional>
#include <mitsuba/core/aabb.h>
#include <mitsuba/mitsuba.h>
#include <optional>
MTS_NAMESPACE_BEGIN

struct MajorantInfo {
  Float majorantDensity;
};

struct SegmentData {
  Float    segDistance;
  uint32_t dataAccessIdx;
};

class Tracker {
public:
  Tracker() = default;

  virtual std::optional<SegmentData> nextSeg() = 0;

  virtual void marchForward(Float step) = 0;

  virtual void marchNext() = 0;

  Float t;
};

class PartitionGrid {
public:
  PartitionGrid() = default;

  PartitionGrid(const AABB &rootBound) : m_rootBound(rootBound){};

  virtual void configure(const std::function<Point3(Point)> &,
                         const std::function<Float(Point3i)> &) = 0;

  virtual void accessData(MajorantInfo *info, uint32_t dataAccessIdx) const = 0;

  virtual std::unique_ptr<Tracker> getTracker(const Ray &ray,
                                              Float      tmax) const = 0;

protected:
  AABB m_rootBound;
};

MTS_NAMESPACE_END