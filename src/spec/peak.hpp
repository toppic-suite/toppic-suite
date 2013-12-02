#ifndef PROT_PEAK_HPP_
#define PROT_PEAK_HPP_

#include <memory>

#include "base/mass_constant.hpp"

namespace prot {

class Peak {
 public:
  Peak(double position, double intensity) {
    position_ = position;
    intensity_ = intensity;
  }

  double getIntensity() {return intensity_;}

  double getPosition() {return position_;}

  void setIntensity(double intensity) {
    intensity_ = intensity;
  }

  void setPosition(double position) {
    position_ = position;
  }

 private:
  double position_;
  double intensity_;
};

inline double compPeakMass(double mono_mz, int charge) {
  return mono_mz * charge - charge * MassConstant::getProtonMass();
}

inline double compMonoMz(double mono_mass, int charge) {
  return mono_mass / charge + MassConstant::getProtonMass();
}

}
#endif
