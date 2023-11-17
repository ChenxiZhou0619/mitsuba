#include "colorspace.h"

namespace spectral {

// clang-format off
const static mitsuba::Matrix3x3 XYZFromsRGB{
    0.4124, 0.3576, 0.1805, 
    0.2126, 0.7152, 0.0722, 
    0.0193, 0.1192, 0.9505
};

const static mitsuba::Matrix3x3 sRGBFromXYZ{
    3.2406, -1.5372, -0.4986,
    -0.9689, 1.8758, 0.0415,
    0.0557, -0.2040, 1.0570
};
// clang-format on

void ColorSpace::Init() {
  ColorSpace::sRGB = std::make_unique<ColorSpace>(sRGBFromXYZ, XYZFromsRGB);
}
} // namespace spectral