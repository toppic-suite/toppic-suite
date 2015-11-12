#ifndef PROT_SPEC_THEO_PEAK_FACTORY_HPP_
#define PROT_SPEC_THEO_PEAK_FACTORY_HPP_

#include "base/activation.hpp"
#include "spec/theo_peak.hpp"

namespace prot {

class TheoPeakFactory {
 public:
  /* called by diagonal.cpp */
  TheoPeakPtrVec geneTheoPeak(BpSpecPtr bp_spec_ptr, ActivationPtr activation_ptr,
                              NeutralLossPtr neutral_loss_ptr,
                              double n_term_shift,double c_term_shift,
                              int bgn,int end,double min_mass, double max_mass);

  TheoPeakPtrVec geneProteoformTheoPeak(ProteoformPtr proteoform_ptr, 
                                        ActivationPtr activation_ptr,
                                        double min_mass);
};

} /* namespace prot */

#endif 
