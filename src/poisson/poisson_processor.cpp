// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "base/logger.hpp"
#include "base/file_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "poisson/poisson_processor.hpp"

//#include <sys/time.h>

namespace prot {

PoissonProcessor::PoissonProcessor(PoissonMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
}

void PoissonProcessor::init() {
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

  comp_ptr_ = PoissonCompPValuePtr(
      new PoissonCompPValue(raw_forms, prot_mod_forms, mng_ptr_));
  LOG_DEBUG("comp pvalue array initialized");
 

  std::string input_file_name = basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->input_file_ext_;
  prsms_ = readPrsm(input_file_name, raw_forms);
  LOG_DEBUG("prsms loaded");
}


/* compute E-value. */

void PoissonProcessor::process() {
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int spectrum_num = countSpNum(spectrum_file_name);
  MsAlignReader reader (spectrum_file_name);
  std::string output_file_name = basename(spectrum_file_name)
       + "." + mng_ptr_->output_file_ext_;
  PrsmWriter writer(output_file_name);
  int cnt = 0;

  /*
  struct timeval start_time; 
  struct timeval end_time; 
  gettimeofday(&start_time, NULL);
  */
  DeconvMsPtr ms_ptr = reader.getNextMs();
  while (ms_ptr.get() != nullptr) {
    cnt++;
    processOneSpectrum(ms_ptr, writer);
    ms_ptr = reader.getNextMs();
    if (ms_ptr.get() != nullptr) {
      std::cout << std::flush << "Poisson computation is processing " << cnt << " of " 
          << spectrum_num << " spectra.\r";
    }
    
    /*
    gettimeofday(&end_time, NULL); 
    float duration = end_time.tv_sec - start_time.tv_sec;
    LOG_DEBUG("Duration: " << duration << " seconds.");
    */
  }
  reader.close();

  //because the prsm_writer ~PrsmWriter changed and the fileclosing is an independant function
  writer.close();
  std::cout << std::endl << "Poisson computation finished." << std::endl;
}

void PoissonProcessor::processOneSpectrum(DeconvMsPtr ms_ptr, PrsmWriter &writer) {
  SpectrumSetPtr spec_set_ptr 
      = getSpectrumSet(ms_ptr, 0, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
  if (spec_set_ptr.get() != nullptr) {
    PrsmPtrVec sele_prsms;
    filterPrsms(prsms_, ms_ptr->getHeaderPtr(), sele_prsms);
    ExtendMsPtr extend_ms_ptr = spec_set_ptr->getMsThreePtr();
    comp_ptr_->setPValueArray(extend_ms_ptr, sele_prsms);
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
