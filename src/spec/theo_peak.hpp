
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

inline bool theoPeakUp(const TheoPeakPtr &a, const TheoPeakPtr &b){
  return a->getPosition() < b->getPosition();
}

std::vector<double> getTheoMassVec (const TheoPeakPtrVec &theo_peaks);

/* called by diagonal.cpp */
TheoPeakPtrVec getTheoPeak(BpSpecPtr bp_spec_ptr, ActivationPtr activation_ptr,
                           NeutralLossPtr neutral_loss_ptr,
                           double n_term_shift,double c_term_shift,
                           int bgn,int end,double min_mass, double max_mass);

TheoPeakPtrVec getProteoformTheoPeak(ProteoformPtr proteoform_ptr, 
                                     ActivationPtr activation_ptr,
                                     double min_mass);
} /* namespace prot */

#endif /* THEO_PEAK_HPP_ */
