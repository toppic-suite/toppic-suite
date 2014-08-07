
#ifndef THEO_PEAK_HPP_
#define THEO_PEAK_HPP_

#include "base/ion.hpp"
#include "base/bp_spec.hpp"
#include "base/activation.hpp"
#include "base/proteoform.hpp"
#include "spec/peak.hpp"

namespace prot {

class TheoPeak : public Peak {
 public:
  TheoPeak(IonPtr ion_ptr,double unmode_mass,double shift);

  IonPtr getIonPtr() {return ion_ptr_;}

  double getModMass() {return getPosition();}

  double getShift() {return shift_;}

 private:
  double unmod_mass_;
  double shift_;
  IonPtr ion_ptr_;
};

typedef std::shared_ptr<TheoPeak> TheoPeakPtr;
typedef std::vector<TheoPeakPtr> TheoPeakPtrVec;

/* called by diagonal.cpp */
TheoPeakPtrVec getTheoPeak(BpSpecPtr bp_spec_ptr, ActivationPtr activation_ptr,
                           NeutralLossPtr neutral_loss_ptr,
                           double n_term_shift,double c_term_shift,
                           int bgn,int end,double min_mass, double max_mass);

TheoPeakPtrVec getProteoformTheoPeak(ProteoformPtr proteoform_ptr, 
                                     ActivationPtr activation_ptr,
                                     double min_mass);

inline bool theoPeakUp(const TheoPeakPtr &a, const TheoPeakPtr &b){
  return a->getPosition() < b->getPosition();
}

/* use inline to speedup */
inline std::vector<double> getTheoMassVec (const TheoPeakPtrVec &theo_peak_list) {
  std::vector<double> masses;
  for (size_t i = 0; i < theo_peak_list.size(); i++) {
    masses.push_back(theo_peak_list[i]->getModMass());
  }
  return masses;
}

} /* namespace prot */

#endif /* THEO_PEAK_HPP_ */
