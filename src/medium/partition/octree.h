/**
 * Octree density partition
 */
#pragma once

#include "partition.h"
#include "tbb/tbb.h"
#include <queue>
MTS_NAMESPACE_BEGIN

struct OcNode {
  bool isLeaf;

  //! If isleaf, it ref to m_majData
  //! If !isLeaf, it ref to the first of sub 8 childrens
  uint32_t childIdx;

  uint32_t parentIdx;
  uint32_t selfIdx;

  Point  nodeCenter;
  Vector nodeSize;

  uint32_t depth;
};

class OctreeGridTracker : public Tracker {
public:
  OctreeGridTracker() = default;

  OctreeGridTracker(Point p, Vector direction, Float t, Float tmax,
                    const OcNode *root);

  virtual std::optional<SegmentData> nextSeg() override;

  virtual void marchForward(Float step) override;

  virtual void marchNext() override;

protected:
  Float cloestCross(Point position) {
    Point pmax = curNode->nodeCenter;

    pmax[0] += ((towards & 0b0100) ? .5f : -.5f) * curNode->nodeSize[0];
    pmax[1] += ((towards & 0b0010) ? .5f : -.5f) * curNode->nodeSize[1];
    pmax[2] += ((towards & 0b0001) ? .5f : -.5f) * curNode->nodeSize[2];

    Float xt = (pmax[0] - position[0]) / direction[0];
    Float yt = (pmax[1] - position[1]) / direction[1];
    Float zt = (pmax[2] - position[2]) / direction[2];

    xt = xt < .0f ? INFINITY : xt;
    yt = yt < .0f ? INFINITY : yt;
    zt = zt < .0f ? INFINITY : zt;

    Float tNear;
    if (xt < yt && xt < zt) {
      stepAxis = 0;
      tNear    = xt;
    } else if (yt < zt) {
      stepAxis = 1;
      tNear    = yt;
    } else {
      stepAxis = 2;
      tNear    = zt;
    }
    stepDistance = tNear;

    return tNear;
  }

private:
  bool terminate = false;

  Point         p;         // current position
  Vector        direction; // march direction
  const OcNode *curNode = nullptr;
  const OcNode *root    = nullptr;

  uint8_t towards;           // along x, y, z axis
  uint8_t stepAxis     = -1; // x = 0, y = 1, z = 2
  Float   stepDistance = .0f;
  Float   tmax;
};

class OctreeGrid : public PartitionGrid {
public:
  OctreeGrid() = delete;

  OctreeGrid(const AABB &rootBound, uint32_t maxDepth, uint32_t minDepth = 4)
      : PartitionGrid(rootBound), m_maxDepth(maxDepth), m_minDepth(minDepth) {
    //
  }

  virtual void
  configure(const std::function<Point3(Point)>  &w2i,
            const std::function<Float(Point3i)> &accessor) override;

  virtual void accessData(MajorantInfo *info,
                          uint32_t      dataAccessIdx) const override;

  virtual std::unique_ptr<Tracker> getTracker(const Ray &ray,
                                              Float      tmax) const override;

private:
  uint32_t m_maxDepth;
  uint32_t m_minDepth;

  std::vector<OcNode> m_octree;
  std::vector<Float>  m_majData;
};

inline OctreeGridTracker::OctreeGridTracker(Point p, Vector direction, Float t,
                                            Float tmax, const OcNode *root) {
  this->p         = p;
  this->t         = t;
  this->root      = root;
  this->direction = direction;
  this->tmax      = tmax;

  // initialize towards
  towards = 0;
  towards |= (direction[0] >= .0f) ? 0b0100 : 0b000;
  towards |= (direction[1] >= .0f) ? 0b0010 : 0b000;
  towards |= (direction[2] >= .0f) ? 0b0001 : 0b000;

  // locate the leaf node contains p
  const OcNode *node = root;
  while (!node->isLeaf) {
    uint32_t offset = 0;
    offset |= (p[0] > node->nodeCenter[0]) ? 0b0100 : 0b0000;
    offset |= (p[1] > node->nodeCenter[1]) ? 0b0010 : 0b0000;
    offset |= (p[2] > node->nodeCenter[2]) ? 0b0001 : 0b0000;

    offset += node->childIdx;
    node = root + offset;
  }

  this->curNode = node;
}

inline std::optional<SegmentData> OctreeGridTracker::nextSeg() {
  if (!curNode || terminate) return std::nullopt;
  Float       distance = cloestCross(p);
  SegmentData seg{distance, curNode->childIdx};
  return std::make_optional<SegmentData>(seg);
}

inline void OctreeGridTracker::marchForward(Float step) {
  t += step;
  p += direction * step;

  if (t >= tmax) terminate = true;
}

inline void OctreeGridTracker::marchNext() {
  marchForward(stepDistance + Epsilon);

  if (terminate) return;

  // locate the next voxel
  // up first
  const OcNode *node = root + curNode->parentIdx;
  // locate the position of curNode in node
  uint32_t offset = curNode->selfIdx - node->childIdx;

  while (true) {
    bool stayInNode = true;
    // pass positive x
    if (stepAxis == 0) {
      stayInNode = (towards & 0b0100) ^ (offset & 0b0100);
    } else if (stepAxis == 1) {
      stayInNode = (towards & 0b0010) ^ (offset & 0b0010);
    } else if (stepAxis == 2) {
      stayInNode = (towards & 0b0001) ^ (offset & 0b0001);
    } else {
      printf("Shouldn't arrive here\n");
      exit(1);
    }

    // up iteratively
    if (!stayInNode) {
      if (node->selfIdx == 0) {
        terminate = true;
        return;
      }
      offset = node->selfIdx - (root + node->parentIdx)->childIdx;
      node   = root + node->parentIdx;
    } else {
      break;
    }
  }

  // stay in node, invert the bit on step axis
  uint32_t base = node->childIdx;
  if (stepAxis == 0)
    offset ^= 0b0100;
  else if (stepAxis == 1)
    offset ^= 0b0010;
  else if (stepAxis == 2)
    offset ^= 0b0001;

  node = root + base + offset;

  // down until leaf node
  curNode = node;
  while (!curNode->isLeaf) {
    uint32_t offset = 0;
    offset |= (p[0] > curNode->nodeCenter[0]) ? 0b0100 : 0b0000;
    offset |= (p[1] > curNode->nodeCenter[1]) ? 0b0010 : 0b0000;
    offset |= (p[2] > curNode->nodeCenter[2]) ? 0b0001 : 0b0000;

    curNode = root + curNode->childIdx + offset;
  }
}

inline void
OctreeGrid::configure(const std::function<Point3(Point)>  &w2i,
                      const std::function<Float(Point3i)> &accessor) {

  std::vector<uint32_t> *leafNodes = new std::vector<uint32_t>();
  // use FIFO to partition the first m_minDepth levels
  std::queue<uint32_t> *toSplit = new std::queue<uint32_t>();

  OcNode rootNode{false,
                  1,
                  0,
                  0,
                  m_rootBound.getCenter(),
                  m_rootBound.max - m_rootBound.min,
                  0};

  // push the rootNode
  m_octree.emplace_back(rootNode);

  toSplit->push(m_octree[0].selfIdx);

  while (!toSplit->empty()) {
    uint32_t nodeIdx = toSplit->front();
    toSplit->pop();

    m_octree[nodeIdx].childIdx = m_octree.size();

    OcNode node    = m_octree[nodeIdx];
    Vector subSize = node.nodeSize * .5f;

    for (int i = 0; i < 8; ++i) {
      OcNode subnode;
      subnode.parentIdx = node.selfIdx;
      subnode.selfIdx   = node.childIdx + i;
      subnode.nodeSize  = subSize;

      Vector centerOffset;
      centerOffset[0]    = ((i & 0b0100) ? .5f : -.5f) * subSize[0];
      centerOffset[1]    = ((i & 0b0010) ? .5f : -.5f) * subSize[1];
      centerOffset[2]    = ((i & 0b0001) ? .5f : -.5f) * subSize[2];
      subnode.nodeCenter = node.nodeCenter + centerOffset;
      subnode.depth      = node.depth + 1;

      if (subnode.depth < m_minDepth)
        // force split before mindepth
        subnode.isLeaf = false;
      else if (subnode.depth >= m_maxDepth)
        // force leaf at maxdepth
        subnode.isLeaf = true;
      else {
        //  split logic here
        Point pmin = subnode.nodeCenter - subnode.nodeSize * .5f;
        Point pmax = subnode.nodeCenter + subnode.nodeSize * .5f;

        auto nullCollisionRate = [&](Point pa, Point pb) {
          //
          Float x0 = INFINITY, y0 = INFINITY, z0 = INFINITY;
          Float x1 = -INFINITY, y1 = -INFINITY, z1 = -INFINITY;

          for (int j = 0; j < 8; ++j) {
            Point pvtx;
            pvtx[0] = (j & 0b0001) ? pa[0] : pb[0];
            pvtx[1] = (j & 0b0010) ? pa[1] : pb[1];
            pvtx[2] = (j & 0b0100) ? pa[2] : pb[2];

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

          Float majorant   = .0f;
          Float times      = .0f;
          Float densitySum = .0f;

          for (int xx = ix0; xx <= ix1; ++xx)
            for (int yy = iy0; yy <= iy1; ++yy)
              for (int zz = iz0; zz <= iz1; ++zz) {
                Float density = accessor({xx, yy, zz});
                majorant      = std::max(majorant, density);
                times += 1.f;
                densitySum += density;
              }

          // null collision rate
          return 1.f - densitySum / (majorant * times);
        };

        Float totalNullCollisionRates = nullCollisionRate(pmin, pmax) * 8.f;
        Float splitNullCollisionRates = .0f;
        for (int j = 0; j < 8; ++j) {
          Point pa = subnode.nodeCenter;

          Vector v = subnode.nodeSize * .5f;
          v[0] *= (j & 0b0100) ? 1.f : -1.f;
          v[1] *= (j & 0b0010) ? 1.f : -1.f;
          v[2] *= (j & 0b0001) ? 1.f : -1.f;
          Point pb = pa + v;

          splitNullCollisionRates += nullCollisionRate(pa, pb);
        }

        if (totalNullCollisionRates > 0.5 ||
            totalNullCollisionRates > 1.2f * splitNullCollisionRates) {
          subnode.isLeaf = false;
        } else {
          subnode.isLeaf = true;
        }
      }

      m_octree.emplace_back(subnode);

      if (!subnode.isLeaf)
        toSplit->push(subnode.selfIdx);
      else
        leafNodes->push_back(subnode.selfIdx);
    }
  }

  delete toSplit;

  int nLeafs = leafNodes->size();
  m_majData  = std::vector<Float>(nLeafs);

  tbb::parallel_for(
      tbb::blocked_range<int>(0, nLeafs), [&](tbb::blocked_range<int> range) {
        for (int r = range.begin(); r < range.end(); ++r) {
          OcNode &leaf  = m_octree[(*leafNodes)[r]];
          leaf.childIdx = r;

          Point pmin = leaf.nodeCenter - leaf.nodeSize * .5f,
                pmax = leaf.nodeCenter + leaf.nodeSize * .5f;

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

          int ix0 = (int)x0;
          int iy0 = (int)y0;
          int iz0 = (int)z0;

          int ix1 = (int)(x1 + 1);
          int iy1 = (int)(y1 + 1);
          int iz1 = (int)(z1 + 1);

          Float majorant = .0f;
          for (int xx = ix0; xx <= ix1; ++xx)
            for (int yy = iy0; yy <= iy1; ++yy)
              for (int zz = iz0; zz <= iz1; ++zz) {
                majorant = std::max(majorant, accessor({xx, yy, zz}));
              }

          m_majData[r] = majorant;
        }
      });

  delete leafNodes;
}

inline void OctreeGrid::accessData(MajorantInfo *info,
                                   uint32_t      dataAccessIdx) const {
  info->majorantDensity = m_majData[dataAccessIdx];
}

inline std::unique_ptr<Tracker> OctreeGrid::getTracker(const Ray &ray,
                                                       Float      tmax) const {
  Float tnear, tfar;
  bool  overlap = m_rootBound.rayIntersect(ray, tnear, tfar);
  if (!overlap || tfar < 0) return nullptr;

  Float tmin = tnear > 0 ? tnear : .0f;
  tmax       = std::min(tmax, tfar);
  return std::make_unique<OctreeGridTracker>(ray(tmin), ray.d, tmin, tmax,
                                             &m_octree[0]);
}

MTS_NAMESPACE_END