//* indirect incident light field

#include "siren.h"
#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class IILF {
public:
  IILF() = default;

  void fit(std::vector<std::pair<Vector, Spectrum>> datas);

private:
  SIREN siren;
};

MTS_NAMESPACE_END