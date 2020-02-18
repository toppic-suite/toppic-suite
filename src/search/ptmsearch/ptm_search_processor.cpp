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
#include "common/util/file_util.hpp"
#include "common/base/prot_mod.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "seq/fasta_reader.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/spectrum_set.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_xml_writer_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "prsm/prsm.hpp"
#include "search/ptmsearch/ptm_search_processor.hpp"
#include "search/ptmsearch/ptm_search_slow_filter.hpp"

namespace toppic {

void seleTopPrsms(const PrsmPtrVec &all_prsm_ptrs, 
                  PrsmPtrVec &sele_prsm_ptrs, int n_report) {
  int match_size = all_prsm_ptrs.size();
  if(all_prsm_ptrs.size()!=0){
    for(int r=0;r< n_report;r++){
      if(r >= match_size){
        break;
      }
      if(all_prsm_ptrs[r]->getMatchFragNum() > 0){
        sele_prsm_ptrs.push_back(all_prsm_ptrs[r]);
      }
    }
  }
  std::sort(sele_prsm_ptrs.begin(), sele_prsm_ptrs.end(), Prsm::cmpSpectrumIdIncPrecursorIdInc);
}

std::function<void()> geneTask(SpectrumSetPtr spectrum_set_ptr, 
                               SimplePrsmPtrVec ori_simple_prsm_ptrs,
                               PtmSearchMngPtr mng_ptr,
                               CompShiftLowMem comp_shift,
                               SimpleThreadPoolPtr pool_ptr, 
                               PrsmXmlWriterSetPtrVec writer_ptr_vec) {

  return [spectrum_set_ptr, ori_simple_prsm_ptrs, mng_ptr, comp_shift, pool_ptr, writer_ptr_vec]() {
    SimplePrsmPtrVec simple_prsm_ptrs = 
        simple_prsm_util::getUniqueMatches(ori_simple_prsm_ptrs);
    PtmSearchSlowFilterPtr slow_filter_ptr = 
        std::make_shared<PtmSearchSlowFilter>(spectrum_set_ptr, simple_prsm_ptrs,
                                              comp_shift, mng_ptr);
    boost::thread::id thread_id = boost::this_thread::get_id();
    int writer_id = pool_ptr->getId(thread_id);
    PrsmXmlWriterSetPtr writer_ptr = writer_ptr_vec[writer_id];

    for (int s = 2; s <= mng_ptr->align_para_ptr_->n_unknown_shift_; s++) {
      PrsmPtrVec complete_prsm_ptrs = slow_filter_ptr->getPrsms(s-2, ProteoformType::COMPLETE);
      std::sort(complete_prsm_ptrs.begin(), complete_prsm_ptrs.end(), 
                Prsm::cmpMatchFragDecStartPosInc);
      PrsmPtrVec sele_complete_prsm_ptrs;
      seleTopPrsms(complete_prsm_ptrs, sele_complete_prsm_ptrs, mng_ptr->n_report_);
      writer_ptr->getCompleteWriterPtr(s)->writeVector(sele_complete_prsm_ptrs);

      PrsmPtrVec prefix_prsm_ptrs = slow_filter_ptr->getPrsms(s-2, ProteoformType::PREFIX);
      std::sort(prefix_prsm_ptrs.begin(), prefix_prsm_ptrs.end(), 
                Prsm::cmpMatchFragDecStartPosInc);
      PrsmPtrVec sele_prefix_prsm_ptrs;
      seleTopPrsms(prefix_prsm_ptrs, sele_prefix_prsm_ptrs, mng_ptr->n_report_);
      writer_ptr->getPrefixWriterPtr(s)->writeVector(sele_prefix_prsm_ptrs);

      PrsmPtrVec suffix_prsm_ptrs = slow_filter_ptr->getPrsms(s-2, ProteoformType::SUFFIX);
      std::sort(suffix_prsm_ptrs.begin(), suffix_prsm_ptrs.end(), Prsm::cmpMatchFragmentDec);
      PrsmPtrVec sele_suffix_prsm_ptrs;
      seleTopPrsms(suffix_prsm_ptrs, sele_suffix_prsm_ptrs, mng_ptr->n_report_);
      writer_ptr->getSuffixWriterPtr(s)->writeVector(sele_suffix_prsm_ptrs);

      PrsmPtrVec internal_prsm_ptrs = slow_filter_ptr->getPrsms(s-2, ProteoformType::INTERNAL);
      std::sort(internal_prsm_ptrs.begin(), internal_prsm_ptrs.end(), Prsm::cmpMatchFragmentDec);
      PrsmPtrVec sele_internal_prsm_ptrs;
      seleTopPrsms(internal_prsm_ptrs, sele_internal_prsm_ptrs, mng_ptr->n_report_);
      writer_ptr->getInternalWriterPtr(s)->writeVector(sele_internal_prsm_ptrs);
    }
  };
}

PtmSearchProcessor::PtmSearchProcessor(PtmSearchMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  comp_shift_ = CompShiftLowMem();
}


// process ptm search
void PtmSearchProcessor::process(){
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = file_util::basename(sp_file_name)+"."+mng_ptr_->input_file_ext_;
  SimplePrsmReader simple_prsm_reader(input_file_name);
  SimplePrsmPtr prsm_ptr = simple_prsm_reader.readOnePrsm();

  // init variables
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  FastaIndexReaderPtr reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);
  int spectrum_num = msalign_util::getSpNum (sp_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  sp_para_ptr->prec_error_ = 0;
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();

  std::string output_file_name = file_util::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          sp_para_ptr->getActivationPtr(),
                          sp_para_ptr->getSkipList());

  const int n_unknown_shift = 2;

  PrsmXmlWriterSetPtrVec writer_set_ptr_vec;
  for (int i = 0; i < mng_ptr_->thread_num_; i++) { 
    std::string writer_file_name = output_file_name + "_" + str_util::toString(i);
    PrsmXmlWriterSetPtr writer_set_ptr = std::make_shared<PrsmXmlWriterSet>(writer_file_name, n_unknown_shift);
    writer_set_ptr_vec.push_back(writer_set_ptr);
  }
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);

  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;
  //LOG_DEBUG("Start search");
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0])!= nullptr){
    cnt+= group_spec_num;
    if(spec_set_ptr->isValid()){
      int spec_id = spec_set_ptr->getSpectrumId();
      SimplePrsmPtrVec selected_prsm_ptrs;
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        selected_prsm_ptrs.push_back(prsm_ptr);
        prsm_ptr = simple_prsm_reader.readOnePrsm();
      }
      if (selected_prsm_ptrs.size() > 0) {
        while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
          boost::this_thread::sleep(boost::posix_time::milliseconds(100));
        }
        pool_ptr->Enqueue(geneTask(spec_set_ptr, selected_prsm_ptrs,
                                   mng_ptr_, comp_shift_, pool_ptr, writer_set_ptr_vec));
      }
    }
    std::cout << std::flush <<  "Multiple PTM search - processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";
  }
  pool_ptr->ShutDown();
  sp_reader.close();
  simple_prsm_reader.close();
  for (int i = 0; i < mng_ptr_->thread_num_; i++) { 
    writer_set_ptr_vec[i]->close();
  }
  LOG_DEBUG("Search completed");
  std::cout << std::endl;

  // Combine results
  int prsm_top_num = mng_ptr_->thread_num_ * mng_ptr_->n_report_;
  for (int s = 2; s <= n_unknown_shift; s++) {
    std::string end_str = "_" + str_util::toString(s);
    // Complete prsms
    std::string complete_output_ext = mng_ptr_->output_file_ext_ + "_" 
        + ProteoformType::COMPLETE->getName() + end_str;
    std::vector<std::string> complete_input_exts;
    for (int t = 0; t < mng_ptr_->thread_num_; t++) {
      std::string input_ext = mng_ptr_->output_file_ext_ + "_" 
          + str_util::toString(t) + "_" + ProteoformType::COMPLETE->getName() + end_str;
      complete_input_exts.push_back(input_ext);
    }
    PrsmStrMergePtr merge_ptr
        = std::make_shared<PrsmStrMerge>(sp_file_name, complete_input_exts, 
                                         complete_output_ext, prsm_top_num);
    merge_ptr->process();
    merge_ptr = nullptr;

    // Prefix prsms
    std::string prefix_output_ext = mng_ptr_->output_file_ext_ + "_" 
        + ProteoformType::PREFIX->getName() + end_str;
    std::vector<std::string> prefix_input_exts;
    for (int t = 0; t < mng_ptr_->thread_num_; t++) {
      std::string input_ext = mng_ptr_->output_file_ext_ + "_" 
          + str_util::toString(t) + "_" + ProteoformType::PREFIX->getName() + end_str;
      prefix_input_exts.push_back(input_ext);
    }
    merge_ptr
        = std::make_shared<PrsmStrMerge>(sp_file_name, prefix_input_exts, 
                                           prefix_output_ext, prsm_top_num);
    merge_ptr->process();
    merge_ptr = nullptr;

    // Suffix prsms
    std::string suffix_output_ext = mng_ptr_->output_file_ext_ + "_" 
        + ProteoformType::SUFFIX->getName() + end_str;
    std::vector<std::string> suffix_input_exts;
    for (int t = 0; t < mng_ptr_->thread_num_; t++) {
      std::string input_ext = mng_ptr_->output_file_ext_ + "_" 
          + str_util::toString(t) + "_" + ProteoformType::SUFFIX->getName() + end_str;
      suffix_input_exts.push_back(input_ext);
    }
    merge_ptr
        = std::make_shared<PrsmStrMerge>(sp_file_name, suffix_input_exts, 
                                         suffix_output_ext, prsm_top_num);
    merge_ptr->process();
    merge_ptr = nullptr;

    // internal prsms
    std::string internal_output_ext = mng_ptr_->output_file_ext_ + "_" 
        + ProteoformType::INTERNAL->getName() + end_str;
    std::vector<std::string> internal_input_exts;
    for (int t = 0; t < mng_ptr_->thread_num_; t++) {
      std::string input_ext = mng_ptr_->output_file_ext_ + "_" 
          + str_util::toString(t) + "_" + ProteoformType::INTERNAL->getName() + end_str;
      internal_input_exts.push_back(input_ext);
    }
    merge_ptr
        = std::make_shared<PrsmStrMerge>(sp_file_name, internal_input_exts, 
                                         internal_output_ext, prsm_top_num);
    merge_ptr->process();
    merge_ptr = nullptr;

    //remove temporary files
    for (int t = 0; t < mng_ptr_->thread_num_; t++) {
      file_util::cleanTempFiles(sp_file_name, complete_input_exts[t]);
      file_util::cleanTempFiles(sp_file_name, prefix_input_exts[t]);
      file_util::cleanTempFiles(sp_file_name, suffix_input_exts[t]);
      file_util::cleanTempFiles(sp_file_name, internal_input_exts[t]);
    }
  }

}

} /* namespace toppic */
