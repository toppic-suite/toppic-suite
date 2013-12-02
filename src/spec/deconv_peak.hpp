#ifndef PROT_DECONV_PEAK_HPP_
#define PROT_DECONV_PEAK_HPP_

#include "spec/peak.hpp"

namespace prot {

class DeconvPeak : public Peak {
 public:
  DeconvPeak (int id, double mono_mass, double intensity, int charge);

  int getCharge() {return charge_;}

  int getId() {return id_;}

  double getMonoMass() {return getPosition();}

  double getMonoMz() {return compMonoMz(getMonoMass(), charge_);}

  double getScore() {return score_;}

  void setId(int id) {id_ = id;}

 private:
  int id_;
  int charge_;
  double score_ = 1.0;
};

typedef std::shared_ptr<DeconvPeak> DeconvPeakPtr;

}
#endif
