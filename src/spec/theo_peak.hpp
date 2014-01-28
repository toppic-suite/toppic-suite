/*
 * theo_peak.hpp
 *
 *  Created on: Nov 28, 2013
 *      Author: xunlikun
 */

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
    TheoPeak(IonPtr ion,double unmode_mass,double shift);
    IonPtr getIonPtr(){return ion_;}
    double getModMass(){return getPosition();}
    double getShift(){return shift_;}
private:
    double unmod_mass_;
    double shift_;

    IonPtr ion_;
};
typedef std::shared_ptr<TheoPeak> TheoPeakPtr;
typedef std::vector<TheoPeakPtr> TheoPeakPtrVec;

inline bool theo_peak_up(const TheoPeakPtr p,TheoPeakPtr n){
  return p->getPosition() < n->getPosition();
}

inline bool theo_peak_down(const TheoPeakPtr p,TheoPeakPtr n){
  return p->getPosition() > n->getPosition();
}

inline void getTheoMassVec (TheoPeakPtrVec &theo_peaks,
                            std::vector<double> &masses) {
  for (unsigned int i = 0; i < theo_peaks.size(); i++) {
    masses.push_back(theo_peaks[i]->getModMass());
  }
}

TheoPeakPtrVec getTheoPeak(BpSpecPtr pep,ActivationPtr type,
                           double n_term_shift,double c_term_shift,
                           int bgn,int end,double min_mass);

TheoPeakPtrVec getProteoformTheoPeak(ProteoformPtr proteoform_ptr, 
                                     ActivationPtr activation_ptr,
                                     double min_mass);
} /* namespace prot */

#endif /* THEO_PEAK_HPP_ */
