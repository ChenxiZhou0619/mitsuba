// clang-format off
#pragma once
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/kdtree.h>
#include "sphericalharmonics.h"
// clang-format on

MTS_NAMESPACE_BEGIN

struct RadianceCacheRecord {
  // store a spherical function
  SphericalHarmonic r, g, b;

  RadianceCacheRecord() : r(6), g(6), b(6) {}
};

struct RadianceRecord {
  Point    position;
  Vector   direction;
  Spectrum radiance;
};

using RadianceNode   = SimpleKDNode<Point, RadianceCacheRecord>;
using RadianceKDTree = PointKDTree<RadianceNode>;

MTS_NAMESPACE_END