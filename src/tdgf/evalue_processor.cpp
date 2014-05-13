#include "base/logger.hpp"
#include "base/file_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "tdgf/evalue_processor.hpp"

namespace prot {

EValueProcessor::EValueProcessor(TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
}

void EValueProcessor::init() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(prsm_para_ptr->getSearchDbFileName(), 
                              prsm_para_ptr->getFixModResiduePtrVec());

  ProteoformPtrVec prot_mod_forms 
      = generateProtModProteoform(raw_forms, prsm_para_ptr->getAllowProtModPtrVec());

  LOG_DEBUG("protein data set loaded");

  ResFreqPtrVec residue_freqs 
      = compResidueFreq(prsm_para_ptr->getFixModResiduePtrVec(), raw_forms); 
  LOG_DEBUG("residue frequency initialized");

  comp_pvalue_ptr_ = CompPValueArrayPtr(
      new CompPValueArray(raw_forms, prot_mod_forms, residue_freqs, mng_ptr_));
  LOG_DEBUG("comp pvalue array initialized");

  std::string input_file_name = basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->input_file_ext_;
  //std::cout<<input_file_name<<std::endl;
  prsms_ = readPrsm(input_file_name, raw_forms);
  LOG_DEBUG("prsms loaded");
}


/* compute E-value. Separate: compute E-value separately or not */

void EValueProcessor::process(bool is_separate) {
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int spectrum_num = countSpNum(spectrum_file_name);
  MsAlignReader reader (spectrum_file_name);
  std::string output_file_name = basename(spectrum_file_name)
       + "." + mng_ptr_->output_file_ext_;
  PrsmWriter writer(output_file_name);
  int cnt = 0;
  DeconvMsPtr ms_ptr = reader.getNextMs();
  while (ms_ptr.get() != nullptr) {
    cnt++;
    processOneSpectrum(ms_ptr, is_separate, writer);
    ms_ptr = reader.getNextMs();
    if (ms_ptr.get() != nullptr) {
      std::cout << std::flush << "E-value computation is processing " << cnt << " of " 
          << spectrum_num << " spectra.\r";
    }
  }
  reader.close();

  //because the prsm_writer ~PrsmWriter changed and the fileclosing is an independant function
  writer.close();
  std::cout << std::endl << "E-value computation finished." << std::endl;
}

void EValueProcessor::processOneSpectrum(DeconvMsPtr ms_ptr, bool is_separate,
                                         PrsmWriter &writer) {
  SpectrumSetPtr spec_set_ptr 
      = getSpectrumSet(ms_ptr, 0, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
  if (spec_set_ptr.get() != nullptr) {
    PrsmPtrVec sele_prsms;
    filterPrsms(prsms_, ms_ptr->getHeaderPtr(), sele_prsms);
    //PrsmUtil.processPrsm(selectedPrsms, deconvSp, seqs);
    
    if (is_separate) {
      for (unsigned i = 0; i < sele_prsms.size(); i++) {
        comp_pvalue_ptr_->setPValue(ms_ptr, sele_prsms[i]);
      }
    } 
    else {
      comp_pvalue_ptr_->setPValueArray(spec_set_ptr->getSpSix(), 
                                       sele_prsms);
    }
    // if matched peak number is too small or E-value is 0, replace it
    // with a max evalue.
    for (unsigned i = 0; i < sele_prsms.size(); i++) {
      if (sele_prsms[i]->getMatchFragNum() < mng_ptr_->comp_evalue_min_match_frag_num_
          || sele_prsms[i]->getEValue() == 0.0) {
        LOG_WARN("Invalid e value!");
        sele_prsms[i]->setProbPtr(getMaxEvaluePtr());
      }
    }
    

    LOG_DEBUG("start sort");
    std::sort(sele_prsms.begin(), sele_prsms.end(), prsmEValueUp);
    LOG_DEBUG("start writing");
    writer.writeVector(sele_prsms);
    LOG_DEBUG("writing complete");
  }
}

}
