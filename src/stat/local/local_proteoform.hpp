//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_STAT_LOCAL_LOCAL_PROTEOFORM_HPP_
#define TOPPIC_STAT_LOCAL_LOCAL_PROTEOFORM_HPP_

#include "seq/mass_shift.hpp"
#include "seq/fasta_seq.hpp"
#include "seq/proteoform.hpp"
#include "stat/local/local_mng.hpp"

namespace toppic {

namespace local_proteoform {

ProteoformPtrVec getCandidateForm(const FastaSeqPtr &seq_ptr, int ori_start_pos, 
                                  int form_start_pos, int form_end_pos, 
                                  MassShiftPtrVec & exp_shift_ptr_vec, 
                                  const LocalMngPtr &mng_ptr); 

ProteoformPtrVec getAllCandidateForms(const ProteoformPtr &ori_form_ptr, 
                                      const LocalMngPtr &mng_ptr); 


ProteoformPtr createProteoformPtr(const ProteoformPtr &base_form_ptr, int match_score, 
                                  std::vector<double> scr_vec, const PtmPtr &ptm_ptr,
                                  const LocalMngPtr &mng_ptr);

ProteoformPtr createProteoformPtr(const ProteoformPtr &base_form_ptr, 
                                  int match_score, int break_pos,
                                  std::vector<double> scr_vec_1, 
                                  std::vector<double> scr_vec_2,
                                  const PtmPtr &ptm_ptr_1, const PtmPtr &ptm_ptr_2,
                                  const LocalMngPtr &mng_ptr);

void getNtermTruncRange(const ProteoformPtr &proteoform, const LocalMngPtr &mng_ptr, int & min, int & max);

void getCtermTruncRange(const ProteoformPtr &proteoform, const LocalMngPtr &mng_ptr, int & min, int & max);

}

}  // namespace toppic

#endif
