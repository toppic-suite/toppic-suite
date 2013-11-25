#ifndef PROT_PEAK_H_
#define PROT_PEAK_H_

#include "mass_constant.hpp"

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

double getPeakMass(double mono_mz, int charge) {
  return mono_mz * charge - charge * MassConstant::getProtonMass();
}

double getMassMonoMz(double mono_mass, int charge) {
  return mono_mass / charge + MassConstant::getProtonMass();
}

}
#endif
