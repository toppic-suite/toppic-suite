/*
 * bp_spec.hpp
 *
 *  Created on: Nov 26, 2013
 *      Author: xunlikun
 */

#ifndef PROT_BP_SPEC_HPP_
#define PROT_BP_SPEC_HPP_

#include "base/residue_seq.hpp"
#include "base/break_point.hpp"

namespace prot {

class BpSpec {
 public:
  BpSpec() {};
  BpSpec(ResSeqPtr res_seq_ptr);
  BreakPointPtrVec getBreakPointPtrVec() {return break_point_ptr_vec_;}
  BreakPointPtr getBreakPointPtr(int i) {return break_point_ptr_vec_[i];}

  /*implement the function getExt[B\Y\C\Z_DOT]masses and getNTermMasses and getCTermMass*/
  std::vector<double> getBreakPointMasses(IonTypePtr ion_type_ptr);

  std::vector<double> getPrmMasses();
  /*implement the function getExt[BY\CZ]masses*/
  std::vector<double> getBreakPointMasses(double n_term_shift,double c_term_shift,
                                          double min_mass,IonTypePtr ion_type_ptr_n,
                                          IonTypePtr ion_type_ptr_c);

  /*implement the function getScaledBMass*/
  std::vector<int> getScaledMass(double scale,IonTypePtr ion_type);
  std::vector<int> getScaledPrmMasses(double scale);
  double getResSeqMass() {return seq_mass_;}
  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  double seq_mass_;
  BreakPointPtrVec break_point_ptr_vec_;

  void initBreakPoints(ResSeqPtr req_seq_ptr);
  /*implement the private function addMass() BpSpec*/
  void addBreakPointMass(double mass,double seq_mass,double min_mass,
                         std::vector<double> &mass_vec);
};

typedef std::shared_ptr<BpSpec> BpSpecPtr;
typedef std::vector<BpSpecPtr> BpSpecPtrVec;

/*come from BpSpec Util*/
int getFirstResPos(double n_term_shift,std::vector<double> ext_b_masses);
int getLastResPos(double c_term_shift,std::vector<double> ext_b_masses);

} /* namespace prot */

#endif /* BP_SPEC_HPP_ */
