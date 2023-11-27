//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#ifndef TOPPIC_FILTER_MNG_VAR_PTM_FILTER_MNG_HPP_
#define TOPPIC_FILTER_MNG_VAR_PTM_FILTER_MNG_HPP_

#include <boost/thread.hpp>

#include "para/prsm_para.hpp"

#include "filter/massmatch/mass_match.hpp"

namespace toppic {

class VarPtmFilterMng {
  public:
  VarPtmFilterMng(PrsmParaPtr prsm_para_ptr,
                  const std::string &index_file_para,
                  const std::string &var_mod_file_name, 
                  int var_ptm_num,
                  int thread_num,
                  bool use_approx_spec,
                  const std::string &output_file_ext);

  static std::vector<int> computeShifts(std::vector<int> &round_single_shift_list,
                                        int var_ptm_num); 

  PrsmParaPtr getPrsmPtr(){return prsm_para_ptr_;}

  std::string getIndexFilePara() {return index_file_para_;}

  int getSingleShiftNum() {return single_shift_list_.size();}

  PrsmParaPtr prsm_para_ptr_;

  std::string index_file_para_;

  //Parameters for fast filteration 
  int max_proteoform_mass_ = 100000;

  //Candidate protein number for each spectrum
  unsigned int comp_num_ = 5;
  unsigned int pref_suff_num_ = 20;
  unsigned int internal_num_ = 30;
  int filter_scale_ = 100;

  int comp_threshold_ = MassMatch::getPrecursorMatchScore() * 2;
  int pref_suff_threshold_ = MassMatch::getPrecursorMatchScore() * 2;
  int internal_threshold_ = MassMatch::getPrecursorMatchScore() * 2 + 4;

  int var_ptm_num_ = 3;

  int thread_num_ = 1;

  bool use_approx_spec_= false;

  int n_spec_block_ = 0;

  std::string output_file_ext_;

  boost::mutex mutex_;

  // processed spectral counts for each database block
  std::vector<int> cnts_;
  
  std::vector<std::string> var_ptm_file_vec_{"zero_ptm_term_index", 
    "zero_ptm_diag_index", "zero_ptm_rev_term_index", "zero_ptm_rev_diag_index"};
  
  // variable PTM mass shift list 
  // use integer to avoid errors in double addition
  double round_scale_ = 10000;
  std::vector<int> round_single_shift_list_;
  std::vector<double> single_shift_list_;

  // shift list for single and multiple (up to var_ptm_num) variable PTMs  
  std::vector<int> round_shift_list_;
  std::vector<double> shift_list_;
  std::vector<int> int_shift_list_;
};

typedef std::shared_ptr<VarPtmFilterMng> VarPtmFilterMngPtr;

} 

#endif 
