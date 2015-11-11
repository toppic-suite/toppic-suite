#ifndef PROT_SPEC_THEO_PEAK_HPP_
#define PROT_SPEC_THEO_PEAK_HPP_

#include "base/ion.hpp"
#include "spec/peak.hpp"

namespace prot {

class TheoPeak;
typedef std::shared_ptr<TheoPeak> TheoPeakPtr;

class TheoPeak : public Peak {
 public:
  TheoPeak(IonPtr ion_ptr,double unmode_mass,double shift);

  IonPtr getIonPtr() {return ion_ptr_;}

  double getModMass() {return getPosition();}

  double getShift() {return shift_;}

  static bool cmpPosIncrease(const TheoPeakPtr &a, const TheoPeakPtr &b){
    return a->getPosition() < b->getPosition();
}

 private:
  IonPtr ion_ptr_;
  double unmod_mass_;
  double shift_;
};

typedef std::vector<TheoPeakPtr> TheoPeakPtrVec;

} /* namespace prot */

#endif /* THEO_PEAK_HPP_ */
