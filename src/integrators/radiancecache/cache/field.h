// clang-format off
#pragma once
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/kdtree.h>
// clang-format on

MTS_NAMESPACE_BEGIN

struct RadianceRecord {
  // store a spherical function
};

using RadianceNode   = SimpleKDNode<Point, RadianceRecord>;
using RadianceKDTree = PointKDTree<RadianceNode>;

MTS_NAMESPACE_END