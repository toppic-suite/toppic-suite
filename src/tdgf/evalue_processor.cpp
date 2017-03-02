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
#include "base/web_logger.hpp"
#include "base/file_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/threadpool.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "tdgf/evalue_processor.hpp"


namespace prot {

void EValueProcessor::init() {
  test_num_ptr_ = std::make_shared<CountTestNum>(mng_ptr_);
  LOG_DEBUG("Count test number initialized.");

  fai_ = fai_load(mng_ptr_->prsm_para_ptr_->getSearchDbFileName().c_str());

  ResFreqPtrVec residue_freqs = test_num_ptr_->getResFreqPtrVec();
  if (!mng_ptr_->use_gf_) {
    comp_pvalue_table_ptr_ = std::make_shared<CompPValueLookupTable>(mng_ptr_);
  }

  comp_pvalue_ptr_ = std::make_shared<CompPValueArray>(test_num_ptr_, mng_ptr_);
  LOG_DEBUG("comp pvalue array initialized");

}

EValueProcessor::~EValueProcessor() {
  fai_destroy(fai_);
}

std::function<void()> geneTask(SpectrumSetPtr spec_set_ptr, const PrsmPtrVec & sele_prsm_ptrs,
                               double ppo, bool is_separate,
                               TdgfMngPtr mng_ptr, CountTestNumPtr test_num_ptr,
                               std::shared_ptr<ThreadPool<PrsmXmlWriter>> pool_ptr) {

  return [spec_set_ptr, sele_prsm_ptrs, ppo, is_separate, mng_ptr, test_num_ptr, pool_ptr]() {
    PrsmPtrVec prsm_vec; // copy sele_prsm_ptrs
    for (size_t i = 0; i < sele_prsm_ptrs.size(); i++) {
      prsm_vec.push_back(std::make_shared<Prsm>(*sele_prsm_ptrs[i].get()));
    }
    CompPValueArrayPtr comp_pvalue_ptr = std::make_shared<CompPValueArray>(test_num_ptr, mng_ptr);
    comp_pvalue_ptr->process(spec_set_ptr, prsm_vec, ppo, is_separate);
    boost::thread::id thread_id = boost::this_thread::get_id();
    PrsmXmlWriterPtr writer_ptr = pool_ptr->getWriter(thread_id);
    for (size_t i = 0; i < prsm_vec.size(); i++) {
      writer_ptr->write(prsm_vec[i]);
    }
  };
}

/* compute E-value. Separate: compute E-value separately or not */
void EValueProcessor::process(bool is_separate) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string spectrum_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string output_file_name = FileUtil::basename(spectrum_file_name) + "." + mng_ptr_->output_file_ext_;
  PrsmXmlWriter writer(output_file_name);

  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();
  std::string input_file_name = FileUtil::basename(spectrum_file_name) + "." + mng_ptr_->input_file_ext_;
  PrsmReader prsm_reader(input_file_name);
  LOG_DEBUG("start read prsm");
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);

  //init variables
  int spectrum_num = MsAlignUtil::getSpNum(spectrum_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(spectrum_file_name, group_spec_num, sp_para_ptr->getActivationPtr());

  std::shared_ptr<ThreadPool<PrsmXmlWriter>> pool_ptr 
      = std::make_shared<ThreadPool<PrsmXmlWriter>>(mng_ptr_->thread_num_ , output_file_name);

  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;

  LOG_DEBUG("Start search");
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    cnt += group_spec_num;
    if(spec_set_ptr->isValid()){
      PrsmPtrVec selected_prsm_ptrs;
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_set_ptr->getSpecId()) {
        selected_prsm_ptrs.push_back(prsm_ptr);
        prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
      }
      if (!mng_ptr_->use_gf_) {
        processOneSpectrum(spec_set_ptr, selected_prsm_ptrs, ppo, is_separate, writer);
      } else if (checkPrsms(selected_prsm_ptrs)) {
        while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ + 2) {
          boost::this_thread::sleep(boost::posix_time::milliseconds(100));
        }
        pool_ptr->Enqueue(geneTask(spec_set_ptr, selected_prsm_ptrs, ppo, is_separate,
                                   mng_ptr_, test_num_ptr_, pool_ptr));
      }
    }

    std::cout << std::flush << "E-value computation - processing " << cnt << " of " 
        << spectrum_num << " spectra.\r";

    if (mng_ptr_->use_gf_){
      WebLog::percentLog(cnt, spectrum_num, WebLog::GfEvalueTime());
    } else {
      WebLog::percentLog(cnt, spectrum_num, WebLog::TableEvalueTime());	
    }
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;
  sp_reader.close();
  prsm_reader.close();
  writer.close();

  PrsmXmlWriterPtr all_writer_ptr = std::make_shared<PrsmXmlWriter>(output_file_name);
  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    PrsmReaderPtr all_reader_ptr = std::make_shared<PrsmReader>(output_file_name + "_" + std::to_string(i)); 
    PrsmPtr p = all_reader_ptr->readOnePrsm(seq_reader, fix_mod_ptr_vec);
    while(p != nullptr) {
      all_writer_ptr->write(p); 
      p = all_reader_ptr->readOnePrsm(seq_reader, fix_mod_ptr_vec); 
    }
    all_reader_ptr->close();
  }
  all_writer_ptr->close();
}

bool EValueProcessor::checkPrsms(const PrsmPtrVec &prsm_ptrs) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    ExtremeValuePtr extreme_value_ptr = prsm_ptrs[i]->getExtremeValuePtr();
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

void EValueProcessor::compEvalues(SpectrumSetPtr spec_set_ptr, PrsmPtrVec &sele_prsm_ptrs,
                                  double ppo, bool is_separate) {

  if (!mng_ptr_->use_gf_ 
      && comp_pvalue_table_ptr_->inTable(spec_set_ptr->getDeconvMsPtrVec(), sele_prsm_ptrs)) {
    comp_pvalue_table_ptr_->process(spec_set_ptr->getDeconvMsPtrVec(), sele_prsm_ptrs, ppo);
    LOG_DEBUG("Using table");
  } else {
    comp_pvalue_ptr_->process(spec_set_ptr, sele_prsm_ptrs, ppo, is_separate);
  }

  // if matched peak number is too small or E-value is 0, replace it
  // with a max evalue.
  for (unsigned i = 0; i < sele_prsm_ptrs.size(); i++) {
    if (sele_prsm_ptrs[i]->getMatchFragNum() <= mng_ptr_->comp_evalue_min_match_frag_num_) {
      sele_prsm_ptrs[i]->setExtremeValuePtr(ExtremeValue::getMaxEvaluePtr());
    } else {
      if (sele_prsm_ptrs[i]->getEValue() == 0.0) {
        LOG_WARN("Invalid e value!");
        sele_prsm_ptrs[i]->setExtremeValuePtr(ExtremeValue::getMaxEvaluePtr());
      }
    }
  }
}

void EValueProcessor::processOneSpectrum(SpectrumSetPtr spec_set_ptr,
                                         PrsmPtrVec &sele_prsm_ptrs,
                                         double ppo, bool is_separate,
                                         PrsmXmlWriter &writer) {
  //LOG_DEBUG("sele prsm number " << sele_prsm_ptrs.size());
  if (spec_set_ptr->isValid()) {

    bool need_comp = checkPrsms(sele_prsm_ptrs);
    //LOG_DEBUG("Need computation: " << need_comp );

    if (need_comp) {
      compEvalues(spec_set_ptr, sele_prsm_ptrs, ppo, is_separate);
    }

    //LOG_DEBUG("start sort");
    std::sort(sele_prsm_ptrs.begin(), sele_prsm_ptrs.end(), Prsm::cmpEValueInc);
    writer.writeVector(sele_prsm_ptrs);
    //LOG_DEBUG("writing complete");
  }
}

}
