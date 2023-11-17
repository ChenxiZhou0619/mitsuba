// clang-format off
#pragma once
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/matrix.h>
// clang-format on

namespace spectral {

class XYZ {
public:
  XYZ() {}
};

class RGB {
  RGB() {}
};

class ColorSpace {
public:
  ColorSpace() = delete;

  ColorSpace(const mitsuba::Matrix3x3 &RGBFromXYZ,
             const mitsuba::Matrix3x3 &XYZFromRGB)
      : m_RGBFromXYZ(RGBFromXYZ), m_XYZFromRGB(XYZFromRGB){};

  static void Init();

  XYZ toXYZ(const RGB &rgb) const;

  RGB toRGB(const XYZ &xyz) const;

public:
  static std::unique_ptr<ColorSpace> sRGB;

private:
  mitsuba::Matrix3x3 m_RGBFromXYZ;
  mitsuba::Matrix3x3 m_XYZFromRGB;
};

} // namespace spectral