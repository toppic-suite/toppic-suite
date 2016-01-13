#ifndef PROT_SPEC_PEAK_UTIL_HPP_
#define PROT_SPEC_PEAK_UTIL_HPP_

#include "base/mass_constant.hpp"

namespace prot {

class PeakUtil {
 public:
  static double compPeakMass(double mono_mz, int charge) {
    return mono_mz * charge - charge * MassConstant::getProtonMass();
  }

  static double compMonoMz(double mono_mass, int charge) {
    return mono_mass / charge + MassConstant::getProtonMass();
  }

};


}
#endif
