#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/thread.h>
#include <mitsuba/render/scene.h>

#include <nanovdb/NanoVDB.h>
#include <nanovdb/util/GridHandle.h>
#include <nanovdb/util/IO.h>
#include <nanovdb/util/SampleFromVoxels.h>

#include "partition/octree.h"
#include "partition/partition.h"
#include "partition/uniformgrid.h"

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

    std::string partitionType = props.getString("partition", "uniformGrid");

    if (partitionType == "uniformGrid") {
      //* Fetch the resolution of uniform grid
      int resX = props.getInteger("majResolutionX", 64);
      int resY = props.getInteger("majResolutionY", 64);
      int resZ = props.getInteger("majResolutionZ", 64);

      m_majGrid =
          std::make_unique<UniformGrid>(m_bounds, Vector3i{resX, resY, resZ});
    } else if (partitionType == "octreeGrid") {
      //* Fetch the maxdepth of octree grid

      m_majGrid = std::make_unique<OctreeGrid>(m_bounds, 7, 2);
    } else {
      Log(EError, "Undefined partition type");
      exit(1);
    }
  }

  virtual void configure() override {
    Medium::configure();

    auto w2i =
        std::bind(&NanovdbMedium::WorldToIndex, this, std::placeholders::_1);
    auto accessor = [&](Point3i p) -> Float {
      const auto density_accessor = m_densityGrid->getAccessor();
      p[0] =
          math::clamp(p[0], m_indexInfo.minIndex[0], m_indexInfo.maxIndex[0]);
      p[1] =
          math::clamp(p[1], m_indexInfo.minIndex[1], m_indexInfo.maxIndex[1]);
      p[2] =
          math::clamp(p[2], m_indexInfo.minIndex[2], m_indexInfo.maxIndex[2]);
      return density_accessor.getValue({p[0], p[1], p[2]});
    };
    m_majGrid->configure(w2i, accessor);
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
                   ETrackingType type, bool *terminated,
                   MajorantSamplingRecord *maj_rec) const override {

    //! Notice, different tracking type using different path-sampling
    //! coefficients

    const uint32_t heroChannel = ray.heroChannel;

    Float    sampledOpacityThick  = -math::fastlog(1 - u);
    Float    accumulatedThick     = .0f;
    Spectrum spectralOpacityThick = Spectrum{.0f};

    // The tracker tracks the majorant grid

    auto tracker = m_majGrid->getTracker(ray, tmax);

    if (!tracker) {
      *terminated          = true;
      maj_rec->free_flight = tmax;
      maj_rec->tr_majorant = Spectrum(1.f);
      maj_rec->pdf_flight  = 1.f;
      return;
    }

    while (auto segOpt = tracker->nextSeg()) {
      auto [segDistance, dataAccessIdx] = *segOpt;

      MajorantInfo info;
      m_majGrid->accessData(&info, dataAccessIdx);

      Float    majorantDensity = info.majorantDensity * m_densityScale;
      Spectrum majorantSigmaT  = majorantDensity * m_sigmaT;

      Float pathSamplingCofficient;
      if (type & ETrackingType::ENaiveDeltaTracking) {
        // Just use the extinction cofficient of hero channel
        pathSamplingCofficient = majorantSigmaT[heroChannel];
      } else if (type & ETrackingType::ESpectralTracking) {
        // Use the extinction cofficient among all channels
        pathSamplingCofficient = majorantSigmaT.max();
      } else {
        Log(EError, "Shouldn't arrive here");
        exit(1);
      }

      Float segThick = segDistance * pathSamplingCofficient;

      if (segThick + accumulatedThick > sampledOpacityThick) {
        // Reach the sampled distance
        Float step =
            (sampledOpacityThick - accumulatedThick) / pathSamplingCofficient;

        spectralOpacityThick += step * majorantSigmaT;
        tracker->marchForward(step);

        Point scatterP = ray(tracker->t);
        // Clamp the density
        Float density = std::min(SampleDensity(scatterP), majorantDensity);

        maj_rec->sigma_a = density * m_sigmaA;
        maj_rec->sigma_s = density * m_sigmaS;
        if (type & ETrackingType::ENaiveDeltaTracking) {
          // Just use the extinction cofficient of hero channel
          maj_rec->sigma_maj   = majorantSigmaT;
          maj_rec->sigma_n     = (majorantDensity - density) * m_sigmaT;
          maj_rec->tr_majorant = (-spectralOpacityThick).exp();

        } else if (type & ETrackingType::ESpectralTracking) {
          // Use the extinction cofficient among all channels
          maj_rec->sigma_maj = Spectrum{pathSamplingCofficient};
          maj_rec->sigma_n =
              Spectrum{pathSamplingCofficient} - density * m_sigmaT;
          maj_rec->tr_majorant = Spectrum{std::exp(-sampledOpacityThick)};
        } else {
          Log(EError, "Shouldn't arrive here");
          exit(1);
        }

        maj_rec->pdf_flight =
            std::exp(-sampledOpacityThick) * pathSamplingCofficient;

        maj_rec->free_flight = tracker->t;
        maj_rec->p           = ray(maj_rec->free_flight);
        maj_rec->medium      = this;

        *terminated = false;
        return;
      }

      // Accumulate the current segment and march to next voxel
      accumulatedThick += segDistance * pathSamplingCofficient;
      spectralOpacityThick += segDistance * majorantSigmaT;
      tracker->marchNext();
    }

    // No nextSeg avaliable, i.e. reach the tmax

    if (type & ETrackingType::ENaiveDeltaTracking) {
      maj_rec->tr_majorant = (-spectralOpacityThick).exp();
    } else if (type & ETrackingType::ESpectralTracking) {
      maj_rec->tr_majorant = Spectrum{std::exp(-accumulatedThick)};
    } else {
      Log(EError, "Shouldn't arrive here");
      exit(1);
    }

    maj_rec->pdf_flight  = std::exp(-accumulatedThick);
    maj_rec->free_flight = tmax;
    maj_rec->p           = ray(maj_rec->free_flight);
    *terminated          = true;
  }

  Float SampleDensity(Point p_world) const override {
    using Sampler =
        nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
    auto p_index = WorldToIndex(p_world);
    return Sampler(m_densityGrid->tree())(p_index) * m_densityScale;
  }

  // FIXME : sample by voxel index is not suitable for all grid type
  virtual Float SampleMajDensity(uint32_t dataAccessIdx) const override {
    MajorantInfo info;
    m_majGrid->accessData(&info, dataAccessIdx);
    return info.majorantDensity;
  };

  virtual std::unique_ptr<Tracker> GetTracker(const Ray &ray,
                                              Float      tmax) const override {
    return m_majGrid->getTracker(ray, tmax);
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

  std::unique_ptr<PartitionGrid> m_majGrid;
};

MTS_IMPLEMENT_CLASS_S(NanovdbMedium, false, Medium)
MTS_EXPORT_PLUGIN(NanovdbMedium, "Nanovdb heterogeneous medium")
MTS_NAMESPACE_END