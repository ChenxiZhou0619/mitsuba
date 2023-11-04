#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/thread.h>
#include <mitsuba/render/scene.h>

#include <nanovdb/NanoVDB.h>
#include <nanovdb/util/GridHandle.h>
#include <nanovdb/util/IO.h>
#include <nanovdb/util/SampleFromVoxels.h>

#include "majorantgrid.h"
#include <tbb/tbb.h>
MTS_NAMESPACE_BEGIN

//* .nvdb format heterogeneous media

using namespace nanovdb;
using BufferT = nanovdb::HostBuffer;

class NanovdbMedium : public Medium {
public:
  NanovdbMedium(const Properties &props) : Medium(props) {
    m_densityScale = props.getFloat("densityScale", 1.f); // Density scale
    m_toWorld      = props.getTransform("toWorld");       // Grid transform
    m_toLocal      = m_toWorld.inverse();

    auto filename = Thread::getThread()->getFileResolver()->resolve(
        props.getString("filename"));
    if (!fs::exists(filename)) {
      Log(EError, "Nanovdbfile %s not exists", filename.c_str());
      exit(1);
    }

    //* Read the density grid
    std::string densityGridName = props.getString("densityGridName", "density");
    m_densityHandle = io::readGrid(filename.c_str(), densityGridName);

    //* Read the temperature grid (optional)
    std::string temperatureGridName =
        props.getString("temperatureGridName", "temperature");
    if (io::hasGrid(filename.c_str(), temperatureGridName)) {
      m_temperaturHandle = io::readGrid(filename.c_str(), temperatureGridName);
    }

    m_densityGrid     = m_densityHandle.grid<Float>();
    m_temperatureGrid = m_temperaturHandle.grid<Float>();

    //* Voxel index bounds
    int minx = m_densityGrid->indexBBox().min().x(),
        miny = m_densityGrid->indexBBox().min().y(),
        minz = m_densityGrid->indexBBox().min().z(),
        maxx = m_densityGrid->indexBBox().max().x(),
        maxy = m_densityGrid->indexBBox().max().y(),
        maxz = m_densityGrid->indexBBox().max().z();

    m_indexInfo.minIndex = {minx, miny, minz};
    m_indexInfo.maxIndex = {maxx, maxy, maxz};

    //* Compute the world bounding box (after transform)
    //* i.e. traverse the 8 vertices
    for (int idx = 0; idx < 8; ++idx) {
      Float ixf      = (idx & 0b0001) ? maxx : minx;
      Float iyf      = (idx & 0b0010) ? maxy : miny;
      Float izf      = (idx & 0b0100) ? maxz : minz;
      Point vtxWorld = IndexToWorld(Point{ixf, iyf, izf});

      m_bounds.expandBy(vtxWorld);
    }

    Log(EInfo,
        "The nanovdbmedium %s 's world bounds [%.2f, %.2f, %.2f], [%.2f, %.2f, "
        "%.2f]",
        filename.c_str(), m_bounds.min.x, m_bounds.min.y, m_bounds.min.z,
        m_bounds.max.x, m_bounds.max.y, m_bounds.max.z);

    //* Fetch the information of majorant grid
    int resX = props.getInteger("majResolutionX", 64);
    int resY = props.getInteger("majResolutionY", 64);
    int resZ = props.getInteger("majResolutionZ", 64);

    m_majResolution = {resX, resY, resZ};
  }

  virtual void configure() override {
    Medium::configure();

    // Create majorant grid
    m_majGrid = std::make_unique<MajorantGrid>(m_bounds, m_majResolution);

    // Initialize majorant grid
    tbb::parallel_for(
        tbb::blocked_range<int>(0, m_majGrid->size()),
        [&](tbb::blocked_range<int> r) {
          for (int i = r.begin(); i < r.end(); ++i) {
            auto [x, y, z]    = m_majGrid->offsetToXYZ(i);
            auto [wmin, wmax] = m_majGrid->voxelBound(x, y, z);
            auto [x0, y0, z0] = WorldToIndex(wmin);
            auto [x1, y1, z1] = WorldToIndex(wmax);

            int ix0 = std::max((int)(x0 - 1), m_indexInfo.minIndex[0]);
            int iy0 = std::max((int)(y0 - 1), m_indexInfo.minIndex[1]);
            int iz0 = std::max((int)(z0 - 1), m_indexInfo.minIndex[2]);

            int ix1 = std::min((int)(x1 + 1), m_indexInfo.maxIndex[0]);
            int iy1 = std::min((int)(y1 + 1), m_indexInfo.maxIndex[1]);
            int iz1 = std::min((int)(z1 + 1), m_indexInfo.maxIndex[2]);

            Float max_density      = .0f;
            auto  density_accessor = m_densityGrid->getAccessor();
            for (int xx = ix0; xx <= ix1; ++xx)
              for (int yy = iy0; yy <= iy1; ++yy)
                for (int zz = iz0; zz <= iz1; ++zz) {
                  max_density = std::max(
                      max_density,
                      density_accessor.getValue({xx, yy, zz}) * m_densityScale);
                }
            m_majGrid->at(x, y, z) = max_density;
          }
        });
  }

  NanovdbMedium(Stream *stream, InstanceManager *manager)
      : Medium(stream, manager) {
    Log(EError, "Nanovdbmedium serilization is not support");
    exit(1);
  }

  void serialize(Stream *stream, InstanceManager *manager) const override {
    Medium::serialize(stream, manager);
    Log(EError, "Nanovdbmedium serilization is not support");
    exit(1);
  }

  Spectrum evalTransmittance(const Ray &ray, Sampler *sampler) const override {
    Log(EError, "vpath doesn't support Nanovdbmedium");
    exit(1);
  }

  bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec,
                      Sampler *sampler) const override {
    Log(EError, "vpath doesn't support Nanovdbmedium");
    exit(1);
  }

  void eval(const Ray &ray, MediumSamplingRecord &mRec) const override {
    Log(EError, "vpath doesn't support Nanovdbmedium");
    exit(1);
  }

  bool isHomogeneous() const override { return false; }

  std::string toString() const override { return "Nanovdbmedium"; }

  virtual void
  sampleTrMajorant(const RayDifferential &ray, Float u, Float tmax,
                   bool                   *terminated,
                   MajorantSamplingRecord *maj_rec) const override {
    // Check if the ray intersect m_bounds, if so, reset the ray origin and tmax

    Float    sampledOpacityThick = -math::fastlog(1 - u);
    Spectrum accumulatedThick    = Spectrum{.0f};

    // The tracker tracks the majorant grid
    DDATracker tracker = GetTracker(ray, tmax);
    while (auto segOpt = tracker.nextSeg()) {

      auto [segDistance, segVoxel] = *segOpt;
      Float majorantDensity =
          m_majGrid->at(segVoxel[0], segVoxel[1], segVoxel[2]);
      Spectrum majorantSigmaT = majorantDensity * m_sigmaT;
      Float    segThick       = segDistance * majorantSigmaT[0];

      if (segThick + accumulatedThick[0] > sampledOpacityThick) {
        // Reach the sampled distance
        Float step =
            (sampledOpacityThick - accumulatedThick[0]) / majorantSigmaT[0];
        accumulatedThick += step * majorantSigmaT;
        tracker.marchForward(step);

        Point scatterP = ray(tracker.t);
        // Clamp the density
        Float density = std::min(SampleDensity(scatterP), majorantDensity);

        maj_rec->sigma_maj   = majorantSigmaT;
        maj_rec->sigma_a     = density * m_sigmaA;
        maj_rec->sigma_s     = density * m_sigmaS;
        maj_rec->sigma_n     = (majorantDensity - density) * m_sigmaT;
        maj_rec->tr_majorant = (-accumulatedThick).exp();

        maj_rec->free_flight = tracker.t;
        maj_rec->p           = ray(maj_rec->free_flight);
        maj_rec->medium      = this;

        *terminated = false;
        return;
      }

      // Accumulate the current segment and march to next voxel
      accumulatedThick += segDistance * majorantSigmaT;
      tracker.marchNext();
    }

    // No nextSeg avaliable, i.e. reach the tmax
    maj_rec->tr_majorant = (-accumulatedThick).exp();

    maj_rec->free_flight = tmax;
    maj_rec->p           = ray(maj_rec->free_flight);
    *terminated          = true;
  }

  MTS_DECLARE_CLASS()

protected:
  //* Transform the index coordinate to worldspace
  Point IndexToWorld(Point p) const {
    // From index to world
    auto xyz = m_densityGrid->indexToWorldF(nanovdb::Vec3f{p[0], p[1], p[2]});
    // Apply user defined transform
    Point world = m_toWorld(Point(xyz[0], xyz[1], xyz[2]));
    return world;
  }

  //* Transform the worldspace coordinate to index
  Point WorldToIndex(Point p) const {
    auto [x, y, z] = m_toLocal(p);
    auto xyz       = m_densityGrid->worldToIndexF(nanovdb::Vec3f{x, y, z});
    return Point{xyz[0], xyz[1], xyz[2]};
  }

  DDATracker GetTracker(const RayDifferential &ray, Float tmax_world) const {
    Float tnear, tfar;
    bool  overlap = m_bounds.rayIntersect(ray, tnear, tfar);
    if (!overlap || tfar < 0) return DDATracker();

    Float tmin = tnear > 0 ? tnear : .0f;
    Float tmax = std::min(tmax_world, tfar);

    Point startVoxel = m_majGrid->worldToIndex(ray(tmin));
    return DDATracker(m_majResolution, m_majGrid->m_voxelSizew, startVoxel,
                      ray.d, tmin, tmax);
  }

  Float SampleDensity(Point p_world) const {
    using Sampler =
        nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
    auto p_index = WorldToIndex(p_world);
    return Sampler(m_densityGrid->tree())(p_index) * m_densityScale;
  }

private:
  //* density scale
  Float m_densityScale = 1.f;

  //* Grid geometry info
  Transform m_toWorld, m_toLocal;
  AABB      m_bounds;

  //* Nanovdb data accessor
  GridHandle<BufferT> m_densityHandle;
  GridHandle<BufferT> m_temperaturHandle;
  const FloatGrid    *m_densityGrid     = nullptr;
  const FloatGrid    *m_temperatureGrid = nullptr;

  //* Grid index info
  struct {
    Point3i minIndex, maxIndex;
  } m_indexInfo;

  //* MajorantGrid
  Vector3i                      m_majResolution; // resolution of majGrid
  std::unique_ptr<MajorantGrid> m_majGrid;
};

MTS_IMPLEMENT_CLASS_S(NanovdbMedium, false, Medium)
MTS_EXPORT_PLUGIN(NanovdbMedium, "Nanovdb heterogeneous medium")
MTS_NAMESPACE_END