//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include <climits>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "seq/proteoform.hpp"
#include "seq/fasta_reader.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_xml_writer_util.hpp"
#include "tdgf/evalue_processor.hpp"


namespace toppic {

void EValueProcessor::init() {
  test_num_ptr_ = std::make_shared<CountTestNum>(mng_ptr_);
  LOG_DEBUG("Count test number initialized.");

  ResFreqPtrVec residue_freqs = test_num_ptr_->getResFreqPtrVec();
  if (!mng_ptr_->use_gf_) {
    comp_pvalue_table_ptr_ = std::make_shared<CompPValueLookupTable>(mng_ptr_);
  }

  comp_pvalue_ptr_ = std::make_shared<CompPValueArray>(test_num_ptr_, mng_ptr_);
  LOG_DEBUG("comp pvalue array initialized");

}

std::function<void()> geneTask(SpectrumSetPtr spec_set_ptr, 
                               const PrsmPtrVec & sele_prsm_ptrs,
                               double ppo, bool is_separate,
                               TdgfMngPtr mng_ptr, CountTestNumPtr test_num_ptr,
                               SimpleThreadPoolPtr pool_ptr,
                               PrsmXmlWriterPtrVec writer_ptr_vec) {

  return [spec_set_ptr, sele_prsm_ptrs, ppo, is_separate, mng_ptr, test_num_ptr, pool_ptr, writer_ptr_vec]() {
    PrsmPtrVec prsm_vec; // copy sele_prsm_ptrs
    for (size_t i = 0; i < sele_prsm_ptrs.size(); i++) {
      prsm_vec.push_back(std::make_shared<Prsm>(*sele_prsm_ptrs[i].get()));
    }
    CompPValueArrayPtr comp_pvalue_ptr = std::make_shared<CompPValueArray>(test_num_ptr, mng_ptr);
    comp_pvalue_ptr->process(spec_set_ptr, prsm_vec, ppo, is_separate);

    for (unsigned i = 0; i < prsm_vec.size(); i++) {
      if (std::round(prsm_vec[i]->getMatchFragNum()) <= std::round(mng_ptr->comp_evalue_min_match_frag_num_)) {
        prsm_vec[i]->setExtremeValuePtr(ExtremeValue::getMaxEvaluePtr());
      } else {
        if (prsm_vec[i]->getEValue() == 0.0) {
          LOG_WARN("Invalid e value!");
          prsm_vec[i]->setExtremeValuePtr(ExtremeValue::getMaxEvaluePtr());
        }
      }
    }

    boost::thread::id thread_id = boost::this_thread::get_id();
    int writer_id = pool_ptr->getId(thread_id);
    for (size_t i = 0; i < prsm_vec.size(); i++) {
      writer_ptr_vec[writer_id]->write(prsm_vec[i]);
    }
  };
}

/* compute E-value. Separate: compute E-value separately or not */
void EValueProcessor::process(bool is_separate) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string output_file_name = file_util::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;
  PrsmXmlWriter writer(output_file_name);

  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();
  std::string input_file_name = file_util::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_;
  PrsmReader prsm_reader(input_file_name);
  LOG_DEBUG("start read prsm");
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);

  // init variables
  int spectrum_num = msalign_util::getSpNum(sp_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          sp_para_ptr->getActivationPtr(),
                          sp_para_ptr->getSkipList());

  PrsmXmlWriterPtrVec writer_ptr_vec = 
      prsm_xml_writer_util::geneWriterPtrVec(output_file_name, mng_ptr_->thread_num_);

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);

  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;

  LOG_DEBUG("Start search");
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0])!= nullptr){
    cnt += group_spec_num;
    if(spec_set_ptr->isValid()){
      PrsmPtrVec selected_prsm_ptrs;
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_set_ptr->getSpectrumId()) {
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
                                   mng_ptr_, test_num_ptr_, pool_ptr, writer_ptr_vec));
      }
    }

    std::cout << std::flush << "E-value computation - processing " << cnt << " of " 
        << spectrum_num << " spectra.\r";
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;
  sp_reader.close();
  prsm_reader.close();
  prsm_xml_writer_util::closeWriterPtrVec(writer_ptr_vec);
  writer.close();

  if (mng_ptr_->use_gf_) {
    int prsm_top_num = INT_MAX; 
    std::vector<std::string> input_exts;
    for (int t = 0; t < mng_ptr_->thread_num_; t++) {
      input_exts.push_back(mng_ptr_->output_file_ext_ + "_" + str_util::toString(t));
    }
    PrsmStrMergePtr merge_ptr
        = std::make_shared<PrsmStrMerge>(sp_file_name, input_exts, 
                                         mng_ptr_->output_file_ext_, prsm_top_num);
    merge_ptr->process();
    merge_ptr = nullptr;
  }

  // remove tempory files
  file_util::cleanTempFiles(sp_file_name, mng_ptr_->output_file_ext_ + "_");
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
    LOG_DEBUG("Fragment number " << sele_prsm_ptrs[i]->getMatchFragNum());
    if (std::round(sele_prsm_ptrs[i]->getMatchFragNum()) <= std::round(mng_ptr_->comp_evalue_min_match_frag_num_)) {
      LOG_DEBUG("Set max e value ");
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
    LOG_DEBUG("sele_prsm_ptrs size " << sele_prsm_ptrs.size());
    writer.writeVector(sele_prsm_ptrs);
    //LOG_DEBUG("writing complete");
  }
}

}
