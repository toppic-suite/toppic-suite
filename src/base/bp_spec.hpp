
#ifndef PROT_BP_SPEC_HPP_
#define PROT_BP_SPEC_HPP_

#include "base/residue_seq.hpp"
#include "base/break_point.hpp"

namespace prot {

// break point spectrum
class BpSpec {
 public:
  BpSpec() {};

  BpSpec(const ResSeqPtr &res_seq_ptr);

  const BreakPointPtrVec& getBreakPointPtrVec() {return break_point_ptr_vec_;}

  BreakPointPtr getBreakPointPtr(int i) {return break_point_ptr_vec_[i];}

  /* Get neutral ion masses for a specific ion type */
  std::vector<double> getBreakPointMasses(IonTypePtr ion_type_ptr);

  std::vector<double> getPrmMasses();

  /* Get rounded scaled neutral ion masses */ 
  std::vector<int> getScaledMass(double scale, IonTypePtr ion_type_ptr);

  std::vector<int> getScaledPrmMasses(double scale);

 private:
  BreakPointPtrVec break_point_ptr_vec_;

  void initBreakPoints(const ResSeqPtr &req_seq_ptr);
};

typedef std::shared_ptr<BpSpec> BpSpecPtr;
typedef std::vector<BpSpecPtr> BpSpecPtrVec;

} /* namespace prot */

#endif /* BP_SPEC_HPP_ */
