#ifndef PROT_GRAPH_ALIGN_MNG_HPP_
#define PROT_GRAPH_ALIGN_MNG_HPP_

#include <string>

#include "prsm/prsm_para.hpp"

namespace prot {

class GraphAlignMng {
 public:
  GraphAlignMng(PrsmParaPtr prsm_para_ptr, 
                const std::string &var_mod_file_name,
                int n_unknown_shift, int max_known_mods,
                int proteo_graph_gap,
                int thread_num,
                const std::string &output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr;
    var_mod_file_name_ = var_mod_file_name;
    n_unknown_shift_ = n_unknown_shift;
    max_known_mods_ = max_known_mods;
    proteo_graph_gap_ = proteo_graph_gap;
    thread_num_ = thread_num;
    output_file_ext_ = output_file_ext;
  }

  PrsmParaPtr prsm_para_ptr_;

  std::string var_mod_file_name_;

  // set it to 1 for testing 
  double error_tolerance_ = 0.1;

  double max_ptm_sum_mass_ = 10000.00;

  double min_consistent_dist_ = 1.0;

  double convert_ratio_ = 274.335215;
  //double convert_ratio_ = 1.0;

  int proteo_graph_gap_;

  int getIntTolerance() {return std::ceil(error_tolerance_ * convert_ratio_);}
  int getIntMaxPtmSumMass() {return std::ceil(max_ptm_sum_mass_ * convert_ratio_);}
  int getIntMinConsistentDist() {return std::ceil(min_consistent_dist_ * convert_ratio_);};

  double align_prefix_suffix_shift_thresh_ = 300;

  double refine_prec_step_width_ = 0.005;

  int max_known_mods_ = 10;

  int n_unknown_shift_ =2;

  int prec_error_ = 2;

  int thread_num_ = 1;

  std::string output_file_ext_;
};

typedef std::shared_ptr<GraphAlignMng> GraphAlignMngPtr;

} /* namespace_prot */

#endif
