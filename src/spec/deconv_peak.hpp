#ifndef PROT_DECONV_PEAK_H_
#define PROT_DECONV_PEAK_H_

#include "peak.hpp"

namespace prot {

class DeconvPeak : public Peak {
 public:
  DeconvPeak (int id, double mono_mass, double intensity, int charge);

  int getCharge() {return charge_;}

  int getId() {return id_;}

  double getMonoMass() {return getPosition();}

  double getMonoMz() {return getMassMonoMz(getMonoMass(), charge_);}

  double getScore() {return score_;}

  void setId(int id) {id_ = id;}

 private:
  int id_;
  int charge_;
  double score_;
};

}
#endif
