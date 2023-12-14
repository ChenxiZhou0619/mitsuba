// clang-format off
#pragma once
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>
#include <mutex>
// clang-format on

MTS_NAMESPACE_BEGIN

namespace pathgrid {

class UniformGrid {
public:
  UniformGrid(const AABB &bound, const Vector3i &res) {
    m_min = bound.min;
    m_max = bound.max;
    m_res = res;

    m_vsize[0] = (m_max[0] - m_min[0]) / res[0];
    m_vsize[1] = (m_max[1] - m_min[1]) / res[1];
    m_vsize[2] = (m_max[2] - m_min[2]) / res[2];

    int totalSize = res.x * res.y * res.z;
    m_reservoirs  = std::vector<Reservoir>(totalSize);
  }

  std::tuple<Vector3i, Vector3i> getNeighborVoxels(Point p, int offset) const {
    Vector to_min = p - m_min;

    int ix_min =
        std::clamp<int>((to_min[0] / m_vsize[0]) - offset, 0, m_res[0] - 1);
    int iy_min =
        std::clamp<int>((to_min[1] / m_vsize[1]) - offset, 0, m_res[1] - 1);
    int iz_min =
        std::clamp<int>((to_min[2] / m_vsize[2]) - offset, 0, m_res[2] - 1);

    int ix_max =
        std::clamp<int>((to_min[0] / m_vsize[0]) + offset, 0, m_res[0] - 1);
    int iy_max =
        std::clamp<int>((to_min[1] / m_vsize[1]) + offset, 0, m_res[1] - 1);
    int iz_max =
        std::clamp<int>((to_min[2] / m_vsize[2]) + offset, 0, m_res[2] - 1);

    return {Vector3i{ix_min, iy_min, iz_min}, Vector3i{ix_max, iy_max, iz_max}};
  }

  std::tuple<Float, Point, Spectrum> getShading(Point p) const {
    Vector to_min = p - m_min;

    int x = std::clamp<int>(to_min[0] / m_vsize[0], 0, m_res[0] - 1);
    int y = std::clamp<int>(to_min[1] / m_vsize[1], 0, m_res[1] - 1);
    int z = std::clamp<int>(to_min[2] / m_vsize[2], 0, m_res[2] - 1);

    int offset = x * m_res[1] * m_res[2] + y * m_res[2] + z;

    Float sw = m_reservoirs[offset].shadingW;
    return {sw, m_reservoirs[offset].shadingP, m_reservoirs[offset].L_y};
  }

  Point getVoxelCenter(Vector3i idx) const {
    Float x = (idx[0] + .5f) * m_vsize[0] + m_min[0];
    Float y = (idx[1] + .5f) * m_vsize[1] + m_min[1];
    Float z = (idx[2] + .5f) * m_vsize[2] + m_min[2];

    return Point{x, y, z};
  }

  void addCandidate(Point p, Vector w, Spectrum L, Float phat, Float pdf,
                    Vector3i idx) {
    int offset = idx[0] * m_res[1] * m_res[2] + idx[1] * m_res[2] + idx[2];
    std::lock_guard<std::mutex> lock_guard(m_reservoirs[offset].mtx);
    m_reservoirs[offset].update(p, w, L, phat, pdf);
  }

  void updateShading() {
    Float M_average = 0;
    Float M_maximum = 0;
    Float M_minimum = 1e20;
    for (int i = 0; i < m_reservoirs.size(); ++i) {
      auto &reservoir = m_reservoirs[i];

      reservoir.updateShadingW();
      //   M_average += reservoir.M;
      //   M_maximum = std::max<Float>(reservoir.M, M_maximum);
      //   M_minimum = std::min<Float>(reservoir.M, M_minimum);

      //   if (reservoir.sum_w > 1e3) {
      //     printf("indx i = %d\n", i);
      //   }
      // }

      // printf("M avg=%.2f, max = %.2f, min = %.2f\n",
      //        M_average / m_reservoirs.size(), M_maximum, M_minimum);
    }
  }

private:
  Point    m_min, m_max; ///< bounding of the uniform grid
  Vector3i m_res;        ///< resolution of the uniform grid
  Vector   m_vsize;      ///< size of each voxel

  struct Reservoir {
    //* Previous reservoir
    int      M_y;
    Float    sum_w_y;
    Float    phat_y;
    Point    p_y;
    Vector   w_y;
    Spectrum L_y;

    //* current reservoir
    int      M_x;
    Float    sum_w_x;
    Float    phat_x;
    Point    p_x;
    Vector   w_x;
    Spectrum L_x;

    std::mutex mtx;

    Float shadingW; ///< previous sum_w/(M * phat_y)
    Point shadingP;

    Reservoir() {
      M_y     = 0;
      sum_w_y = .0f;

      M_x     = 0;
      sum_w_x = .0f;

      shadingW = .0f;
    }

    void update(Point p, Vector w, Spectrum L, Float phat, Float pdf) {
      M_x++;
      Float weight_x = phat / pdf;
      sum_w_x += weight_x;

      if ((Float)rand() / RAND_MAX < (weight_x / sum_w_x)) {
        // update
        p_x    = p;
        w_x    = w;
        L_x    = L;
        phat_x = phat;
      }
    }

    void updateShadingW() {
      // Combine the previous reservoir and current reservoir

      // clamp if M_y is too high
      if (M_y >= 10 * M_x) {
        int new_M_y = 10 * M_x;
        sum_w_y     = sum_w_y / M_y * new_M_y;
        M_y         = new_M_y + M_x;
      } else {
        M_y += M_x;
      }

      // resample between reservoirs
      sum_w_y += sum_w_x;
      if (((Float)rand() / RAND_MAX) < (sum_w_x / sum_w_y)) {
        p_y    = p_x;
        w_y    = w_x;
        L_y    = L_x;
        phat_y = phat_x;
      }

      // compute shaingw
      shadingW = (M_y == 0) ? .0f : (sum_w_y / (M_y * phat_y));
      shadingP = p_y;

      // clean current reservoir
      M_x     = 0;
      sum_w_x = .0f;
    }
  };

  std::vector<Reservoir> m_reservoirs;
};

} // namespace pathgrid

MTS_NAMESPACE_END
