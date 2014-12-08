//#include <sys/time.h>

#include "base/logger.hpp"
#include "base/file_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "tdgf/evalue_processor.hpp"


namespace prot {

EValueProcessor::EValueProcessor(TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
}

void EValueProcessor::init() {
  CountTestNumPtr test_num_ptr = CountTestNumPtr(new CountTestNum(mng_ptr_));
  LOG_DEBUG("Count test number initialized.");

  fai_ = fai_load(mng_ptr_->prsm_para_ptr_->getSearchDbFileName().c_str());

  ResFreqPtrVec residue_freqs = test_num_ptr->getResFreqPtrVec();
  if (mng_ptr_->use_table) {
    comp_pvalue_table_ptr_ = CompPValueLookupTablePtr(
        new CompPValueLookupTable(mng_ptr_));
  }

  comp_pvalue_ptr_ = CompPValueArrayPtr(new CompPValueArray(test_num_ptr, mng_ptr_));
  LOG_DEBUG("comp pvalue array initialized");

}

EValueProcessor::~EValueProcessor() {
  fai_destroy(fai_);
}


/* compute E-value. Separate: compute E-value separately or not */

void EValueProcessor::process(bool is_separate) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string spectrum_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = basename(spectrum_file_name)+"."+mng_ptr_->input_file_ext_;
  PrsmReader prsm_reader(input_file_name);
  std::string output_file_name = basename(spectrum_file_name) + "." + mng_ptr_->output_file_ext_;
  PrsmWriter writer(output_file_name);

  ResiduePtrVec residue_ptr_vec = prsm_para_ptr->getFixModResiduePtrVec();
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fai_, residue_ptr_vec);

  //init variables
  int spectrum_num = getSpNum (spectrum_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(spectrum_file_name, group_spec_num);
  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;

  LOG_DEBUG("Start search");
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    cnt+= group_spec_num;
    if(spec_set_ptr->isValid()){
      PrsmPtrVec selected_prsm_ptrs;
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_set_ptr->getSpecId()) {
        selected_prsm_ptrs.push_back(prsm_ptr);
        prsm_ptr = prsm_reader.readOnePrsm(fai_, residue_ptr_vec);
      }
      processOneSpectrum(spec_set_ptr, selected_prsm_ptrs, is_separate, writer);
    }

    std::cout << std::flush << "E-value computation is processing " << cnt << " of " 
        << spectrum_num << " spectra.\r";

    if (mng_ptr_->use_table == false){
      WebLog::percent_log(0.26 + (double) cnt / spectrum_num * 0.73);
    } else {
      WebLog::percent_log((0.91 + (double) cnt / spectrum_num * 0.07) / 3.5);	
    }
  }
  std::cout << std::endl;
  sp_reader.close();
  prsm_reader.close();
  writer.close();
}

bool EValueProcessor::checkPrsms(const PrsmPtrVec &prsm_ptrs) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    ExtremeValuePtr extreme_value_ptr = prsm_ptrs[i]->getProbPtr();
    if (extreme_value_ptr != nullptr) {
      double evalue = extreme_value_ptr->getEValue();
      double frag_num = prsm_ptrs[i]->getMatchFragNum();
      if (evalue <= mng_ptr_->computation_evalue_cutoff
          && frag_num >= mng_ptr_->computation_frag_num_cutoff) {
        return false;
      }
    }
  }
  return true;
}

void EValueProcessor::compEvalues(SpectrumSetPtr spec_set_ptr,
                                  bool is_separate, PrsmPtrVec &sele_prsm_ptrs) {
<<<<<<< HEAD

  if (mng_ptr_->use_table && comp_pvalue_table_ptr_->inTable(spec_set_ptr->getDeconvMsPtr(), sele_prsm_ptrs)) {
    comp_pvalue_table_ptr_->process(spec_set_ptr->getDeconvMsPtr(), sele_prsm_ptrs);
    LOG_DEBUG("Using table");
=======
  
  if (mng_ptr_->use_table) {
    comp_pvalue_table_ptr_->process(spec_set_ptr->getDeconvMsPtrVec(), sele_prsm_ptrs);
>>>>>>> Completed tdgf for multiple spectra
  }
  else {
    comp_pvalue_ptr_->process(spec_set_ptr, is_separate, sele_prsm_ptrs);
  }

  // if matched peak number is too small or E-value is 0, replace it
  // with a max evalue.
  for (unsigned i = 0; i < sele_prsm_ptrs.size(); i++) {
    if (sele_prsm_ptrs[i]->getMatchFragNum() <= mng_ptr_->comp_evalue_min_match_frag_num_) {
      sele_prsm_ptrs[i]->setProbPtr(getMaxEvaluePtr());
    }
    else {
      if (sele_prsm_ptrs[i]->getEValue() == 0.0) {
        LOG_WARN("Invalid e value!");
        sele_prsm_ptrs[i]->setProbPtr(getMaxEvaluePtr());
      }
    }
  }
}

<<<<<<< HEAD
void EValueProcessor::processOneSpectrum(DeconvMsPtr ms_ptr,
=======
void EValueProcessor::processOneSpectrum(SpectrumSetPtr spec_set_ptr,
>>>>>>> Completed tdgf for multiple spectra
                                         PrsmPtrVec &sele_prsm_ptrs,
                                         bool is_separate,
                                         PrsmWriter &writer) {
  //LOG_DEBUG("sele prsm number " << sele_prsm_ptrs.size());
<<<<<<< HEAD
  SpectrumSetPtr spec_set_ptr
      = getSpectrumSet(ms_ptr, 0, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
  if (spec_set_ptr.get() != nullptr) {
=======
  if (spec_set_ptr->isValid()) {
>>>>>>> Completed tdgf for multiple spectra

    bool need_comp = checkPrsms(sele_prsm_ptrs);
    //LOG_DEBUG("Need computation: " << need_comp );

    if (need_comp) {
      compEvalues(spec_set_ptr, is_separate, sele_prsm_ptrs);
    }

    //LOG_DEBUG("start sort");
    std::sort(sele_prsm_ptrs.begin(), sele_prsm_ptrs.end(), prsmEValueUp);
    writer.writeVector(sele_prsm_ptrs);
    //LOG_DEBUG("writing complete");
  }
}

}
