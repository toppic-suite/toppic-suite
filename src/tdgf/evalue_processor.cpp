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
  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(mng_ptr_->search_db_file_name_,
                              mng_ptr_->fix_mod_residue_list_);

  ProteoformPtrVec prot_mod_forms 
      = generateProtModProteoform(raw_forms, mng_ptr_->allow_prot_mod_list_); 

  LOG_DEBUG("protein data set loaded");

  ResFreqPtrVec residue_freqs 
      = compResidueFreq(mng_ptr_->fix_mod_residue_list_, raw_forms); 
  LOG_DEBUG("residue frequency initialized");

  comp_pvalue_ptr_ = CompPValueArrayPtr(
      new CompPValueArray(raw_forms, prot_mod_forms, residue_freqs, mng_ptr_));
  LOG_DEBUG("comp pvalue array initialized");

  std::string input_file_name = basename(mng_ptr_->spectrum_file_name_)
      + "." + mng_ptr_->input_file_ext_;
  std::cout<<input_file_name<<std::endl;
  prsms_ = readPrsm(input_file_name, raw_forms);
  LOG_DEBUG("prsms loaded");
}


/* compute E-value. Separate: compute E-value separately or not */

void EValueProcessor::process(bool is_separate) {
  int spectrum_num = countSpNum(mng_ptr_->spectrum_file_name_);
  MsAlignReader reader (mng_ptr_->spectrum_file_name_);
  std::string output_file_name = basename(mng_ptr_->spectrum_file_name_)
       + "." + mng_ptr_->output_file_ext_;
  PrSMWriter writer(output_file_name);
  int cnt = 0;
  DeconvMsPtr ms_ptr = reader.getNextMs();
  while (ms_ptr.get() != nullptr) {
    cnt++;
    processOneSpectrum(ms_ptr, is_separate, writer);
    ms_ptr = reader.getNextMs();
    if (ms_ptr.get() != nullptr) {
      std::cout << std::flush << "Evalue computation complete " << cnt << " of " 
          << spectrum_num << " scan number " << ms_ptr->getHeaderPtr()->getScansString() << std::endl;
    }
  }
  reader.close();

  //because the prsm_writer ~PrSMWriter changed and the fileclosing is an independant function
  writer.close();
}

void EValueProcessor::processOneSpectrum(DeconvMsPtr ms_ptr, bool is_separate,
                                         PrSMWriter &writer) {
  SpectrumSetPtr spec_set_ptr 
      = getSpectrumSet(ms_ptr, 0, mng_ptr_->sp_para_ptr_);
  if (spec_set_ptr.get() != nullptr) {
    PrSMPtrVec sele_prsms;
    filterPrsms(prsms_, ms_ptr->getHeaderPtr(), sele_prsms);
    //PrSMUtil.processPrSM(selectedPrsms, deconvSp, seqs);
    
    if (is_separate) {
      for (unsigned i = 0; i < sele_prsms.size(); i++) {
        comp_pvalue_ptr_->setPValue(ms_ptr, sele_prsms[i]);
      }
    } 
    else {
      comp_pvalue_ptr_->setPValueArray(spec_set_ptr->getSpSix(), 
                                       sele_prsms);
    }
    LOG_DEBUG("start sort");
    std::sort(sele_prsms.begin(), sele_prsms.end(), prsmEValueUp);
    LOG_DEBUG("start writing");
    writer.writeVector(sele_prsms);
    LOG_DEBUG("writing complete");
  }
}

}
