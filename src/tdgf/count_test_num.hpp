#ifndef PROT_COUNT_TEST_NUM_HPP_
#define PROT_COUNT_TEST_NUM_HPP_

#include "base/proteoform.hpp"
#include "tdgf/tdgf_mng.hpp"

namespace prot {

class CountTestNum {
 public:
  CountTestNum();

 private:
  static double PREFIX_SUFFIX_ADJUST() {return 0.693;}
  static double INTERNAL_ADJUST() {return 0.508;}
	
	TdgfMngPtr mng_ptr_;

  ProteoformPtrVec raw_forms_;
  ProteoformPtrVec prot_mod_forms_;
	
	double *comp_mass_cnts_;
	double *prec_mass_cnts_;
	double *suff_mass_cnts_;
	double *internal_mass_cnts_;

	double convert_ratio_;
	int max_sp_len_;
	int residue_avg_len_;
	double norm_factor_;
};

}

#endif
