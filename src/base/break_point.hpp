#ifndef PROT_BASE_BREAK_POINT_HPP_
#define PROT_BASE_BREAK_POINT_HPP_

#include "base/ion_type.hpp"

namespace prot {

class BreakPoint {
 public:
  BreakPoint(double prm,double srm);

  double getPrm(){return prm_;}

  double getSrm(){return srm_;}

  double getNTermMass(IonTypePtr ion_type_ptr) {
    return prm_ + ion_type_ptr->getShift();}

  double getCTermMass(IonTypePtr ion_type_ptr) {
    return srm_ + ion_type_ptr->getShift();}

 private:
  double prm_;
  double srm_;
};

typedef std::shared_ptr<BreakPoint> BreakPointPtr;
typedef std::vector<BreakPointPtr> BreakPointPtrVec;

} /* namespace prot */

#endif /* BREAKPOINT_HPP_ */
