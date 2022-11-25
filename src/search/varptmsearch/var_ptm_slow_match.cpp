//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "ms/factory/prm_ms_util.hpp"
#include "search/diag/diag_header_util.hpp"
#include "search/diag/diag_pair_util.hpp"
#include "search/varptmsearch/var_ptm_slow_match.hpp"

namespace toppic {

VarPtmSlowMatch::VarPtmSlowMatch(ProteoformPtr proteo_ptr,
                                 SpectrumSetPtr spectrum_set_ptr,
                                 VarPtmSearchMngPtr mng_ptr) {
  proteo_ptr_ = proteo_ptr;
  deconv_ms_ptr_vec_ = spectrum_set_ptr->getDeconvMsPtrVec();
  ms_six_ptr_vec_ = spectrum_set_ptr->getMsSixPtrVec();
  prec_mono_mass_ = spectrum_set_ptr->getPrecMonoMass();
  mng_ptr_ = mng_ptr;
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  prec_error_tole_ = deconv_ms_ptr_vec_[0]->getMsHeaderPtr()->getPrecErrorTolerance(tole_ptr->getPpo());
  init();
}

inline DiagHeaderPtrVec VarPtmSlowMatch::geneVarPtmNTermShiftHeaders() {
  DiagHeaderPtrVec header_ptrs;

  double seq_mass = proteo_ptr_->getResSeqPtr()->getSeqMass();
  bool exist_n_term = false;
  bool exist_c_term = false;
  // add shift masses
  std::vector<double> shifts = mng_ptr_->shift_list_;
  for (size_t i = 0; i < shifts.size(); i++) {
    double n_shift = shifts[i]; 
    // n_term strict; c_term nostrict; prot n_term no_match; prot c_term no_match
    // pep n_term match; pep c_term no_match
    DiagHeaderPtr header_ptr 
      = std::make_shared<DiagHeader>(n_shift, true, false, false, false, true, false);
    // find Protein N-terminal and C-terminal matches 
    if (n_shift == 0.0) {
      header_ptr->setProtNTermShift(true);
      exist_n_term = true;
    }
    else {
      header_ptr->setProtNTermShift(false);
    }

    double c_shift = prec_mono_mass_ - seq_mass - n_shift;
    if (std::abs(c_shift) <= prec_error_tole_) {
      header_ptr->setProtCTermShift(true);
      exist_c_term = true;
    } 
    else {
      header_ptr->setProtCTermShift(false);
    }
    header_ptr->initHeader(c_shift, proteo_ptr_); 
    header_ptrs.push_back(header_ptr);
  }
  if (!exist_n_term || !exist_c_term) {
    LOG_ERROR("Error: No terminal shifts are found!");
    exit(EXIT_FAILURE);
  }
  return header_ptrs;
}

// initialize var_align
void VarPtmSlowMatch::init() {
  DiagHeaderPtrVec n_term_shift_header_ptrs = geneVarPtmNTermShiftHeaders(); 
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  PrmPeakPtrVec prm_peaks = prm_ms_util::getPrmPeakPtrs(ms_six_ptr_vec_, tole_ptr);
  int group_spec_num = ms_six_ptr_vec_.size();
  DiagonalPtrVec diagonal_ptrs = diag_pair_util::geneDiagonals(n_term_shift_header_ptrs,
                                                               prm_peaks, group_spec_num,
                                                               proteo_ptr_);

  std::vector<double> seq_masses = proteo_ptr_->getBpSpecPtr()->getPrmMasses();
  std::vector<double> ms_masses(prm_peaks.size());
  for (size_t i = 0; i < prm_peaks.size(); i++) {
    ms_masses[i] = prm_peaks[i]->getPosition();
  }
  int bgn = proteo_ptr_->getStartPos();
  int end = proteo_ptr_->getEndPos();
  ResSeqPtr sub_res_seq_ptr = proteo_ptr_->getResSeqPtr()->getSubResidueSeq(bgn, end);
  var_ptm_align_ptr_ = std::make_shared<VarPtmAlign>(ms_masses, seq_masses, 
                                                     sub_res_seq_ptr,  
                                                     diagonal_ptrs, mng_ptr_);
}


PrsmPtr VarPtmSlowMatch::compute() {
  var_ptm_align_ptr_->compute(); 
  return var_ptm_align_ptr_->geneResult(proteo_ptr_, deconv_ms_ptr_vec_,
                                        mng_ptr_->prsm_para_ptr_);
}


}  // namespace toppic
