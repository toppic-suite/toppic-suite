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

#include <string>

#include "base/proteoform.hpp"
#include "base/proteoform_factory.hpp"
#include "base/file_util.hpp"
#include "base/thread_pool.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "zeroptmfilter/zero_ptm_filter_processor.hpp"
#include "zeroptmfilter/mass_zero_ptm_filter.hpp"

namespace prot {

void ZeroPtmFilterProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader ms_reader(prsm_para_ptr->getSpectrumFileName(),
                          group_spec_num,
                          prsm_para_ptr->getSpParaPtr()->getActivationPtr(),
                          prsm_para_ptr->getSpParaPtr()->getSkipList());
  int thread_num = mng_ptr_->thread_num_;
  std::vector<std::ofstream *> output_vec;
  for (int i = 0; i < thread_num; i++) {
    output_vec.push_back(new std::ofstream(sp_file_name + "_" + std::to_string(i)));
  }

  std::vector<std::string> ms_lines = ms_reader.readOneSpectrum();
  int cnt = 0;
  while (ms_lines.size() > 0) {
    cnt = cnt % thread_num;
    for (int i = 0; i < group_spec_num; i++) {
      (*output_vec[cnt]) << std::endl;
      for (size_t k = 0; k < ms_lines.size(); k++) {
        (*output_vec[cnt]) << ms_lines[k] << std::endl;
      }
      (*output_vec[cnt]) << std::endl;
      ms_lines = ms_reader.readOneSpectrum();
    }
    cnt++;
  }

  ms_reader.close();

  for (size_t i = 0; i < output_vec.size(); i++) {
    output_vec[i]->close();
  }

  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    std::cout << "Non PTM filtering - block " << (i + 1) << " out of "
        << db_block_ptr_vec.size() << " started." << std::endl;
    processBlock(db_block_ptr_vec[i]);
    std::cout << "Non PTM filtering - block " << (i + 1) << " finished. " << std::endl;
  }

  std::cout << "Non PTM filtering - combining blocks started." << std::endl;

  int block_num = db_block_ptr_vec.size();

  SimplePrsmStrCombine comp_combine(sp_file_name, mng_ptr_->output_file_ext_ + "_COMPLETE",
                                    block_num, mng_ptr_->output_file_ext_ + "_COMPLETE",
                                    mng_ptr_->comp_num_);
  comp_combine.process();

  SimplePrsmStrCombine pref_combine(sp_file_name, mng_ptr_->output_file_ext_ + "_PREFIX",
                                    block_num, mng_ptr_->output_file_ext_ + "_PREFIX",
                                    mng_ptr_->pref_suff_num_);
  pref_combine.process();

  SimplePrsmStrCombine suff_combine(sp_file_name, mng_ptr_->output_file_ext_ + "_SUFFIX",
                                    block_num, mng_ptr_->output_file_ext_ + "_SUFFIX",
                                    mng_ptr_->pref_suff_num_);
  suff_combine.process();

  SimplePrsmStrCombine internal_combine(sp_file_name, mng_ptr_->output_file_ext_ + "_INTERNAL",
                                        block_num, mng_ptr_->output_file_ext_ + "_INTERNAL",
                                        mng_ptr_->inte_num_);
  internal_combine.process();

  std::cout << "Non PTM filtering - combining blocks finished." << std::endl;
}

std::function<void()> geneTask(const ProteoformPtrVec & raw_forms,
                               const std::string & block_str,
                               ZeroPtmFilterMngPtr mng_ptr, int idx) {
  return[raw_forms, block_str, mng_ptr, idx] () {
    int group_spec_num = mng_ptr->prsm_para_ptr_->getGroupSpecNum();
    MassZeroPtmFilterPtr filter_ptr = std::make_shared<MassZeroPtmFilter>(raw_forms, mng_ptr);
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
    MsAlignReader reader(prsm_para_ptr->getSpectrumFileName() + "_" + std::to_string(idx),
                         group_spec_num,
                         prsm_para_ptr->getSpParaPtr()->getActivationPtr(),
                         prsm_para_ptr->getSpParaPtr()->getSkipList());
    std::string output_file_name = file_util::basename(prsm_para_ptr->getSpectrumFileName())
        + "." + mng_ptr->output_file_ext_;

    SimplePrsmXmlWriter comp_writer(output_file_name + "_COMPLETE_" + block_str + "_" + std::to_string(idx));
    SimplePrsmXmlWriter pref_writer(output_file_name + "_PREFIX_" + block_str + "_" + std::to_string(idx));
    SimplePrsmXmlWriter suff_writer(output_file_name + "_SUFFIX_" + block_str + "_" + std::to_string(idx));
    SimplePrsmXmlWriter internal_writer(output_file_name + "_INTERNAL_" + block_str + "_" + std::to_string(idx));

    std::vector<SpectrumSetPtr> spec_set_vec = reader.getNextSpectrumSet(sp_para_ptr);
    int spectrum_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
    int cnt = 0;
    while (spec_set_vec[0] != nullptr) {
      cnt+= group_spec_num * mng_ptr->thread_num_;
      for (size_t k = 0; k < spec_set_vec.size(); k++) {
        LOG_DEBUG("spec set ptr valid " << spec_set_vec[k]->isValid());
        if (spec_set_vec[k]->isValid()) {
          ExtendMsPtrVec ms_ptr_vec = spec_set_vec[k]->getMsThreePtrVec();
          filter_ptr->computeBestMatch(ms_ptr_vec);
          comp_writer.write(filter_ptr->getCompMatchPtrs());
          pref_writer.write(filter_ptr->getPrefMatchPtrs());
          suff_writer.write(filter_ptr->getSuffMatchPtrs());
          internal_writer.write(filter_ptr->getInternalMatchPtrs());
        }
      }
      if (idx == 0) {
        std::cout << std::flush << "Non PTM filtering - processing " << cnt
            << " of " << spectrum_num << " spectra.\r";
      }
      spec_set_vec = reader.getNextSpectrumSet(sp_para_ptr);
    }
    reader.close();
    comp_writer.close();
    pref_writer.close();
    suff_writer.close();
    internal_writer.close();
  };
}

void ZeroPtmFilterProcessor::processBlock(DbBlockPtr block_ptr) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  int spectrum_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
      + "_" + std::to_string(block_ptr->getBlockIdx());
  ProteoformPtrVec raw_forms
      = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                        prsm_para_ptr->getFixModPtrVec());
  std::string block_str = std::to_string(block_ptr->getBlockIdx());

  std::vector<ThreadPtr> thread_vec;
  for (int i = 1; i < mng_ptr_->thread_num_; i++) {
    ThreadPtr thread_ptr = std::make_shared<boost::thread>(geneTask(raw_forms, block_str, mng_ptr_, i));
    thread_vec.push_back(thread_ptr);
  }

  std::function<void()> task = geneTask(raw_forms, block_str, mng_ptr_, 0);
  task();

  for (size_t i = 0; i < thread_vec.size(); i++) {
    if (thread_vec[i]->joinable()) thread_vec[i]->join();
  }

  std::cout << std::flush << "Non PTM filtering - processing " << spectrum_num
      << " of " << spectrum_num << " spectra." << std::endl;

  SimplePrsmStrCombine comp_combine(sp_file_name,
                                    mng_ptr_->output_file_ext_ + "_COMPLETE_" + block_str,
                                    mng_ptr_->thread_num_,
                                    mng_ptr_->output_file_ext_ + "_COMPLETE_" + block_str,
                                    mng_ptr_->comp_num_);
  comp_combine.process();

  SimplePrsmStrCombine pref_combine(sp_file_name,
                                    mng_ptr_->output_file_ext_ + "_PREFIX_" + block_str,
                                    mng_ptr_->thread_num_,
                                    mng_ptr_->output_file_ext_ + "_PREFIX_" + block_str,
                                    mng_ptr_->pref_suff_num_);
  pref_combine.process();

  SimplePrsmStrCombine suff_combine(sp_file_name,
                                    mng_ptr_->output_file_ext_ + "_SUFFIX_" + block_str,
                                    mng_ptr_->thread_num_,
                                    mng_ptr_->output_file_ext_ + "_SUFFIX_" + block_str,
                                    mng_ptr_->pref_suff_num_);
  suff_combine.process();

  SimplePrsmStrCombine internal_combine(sp_file_name,
                                        mng_ptr_->output_file_ext_ + "_INTERNAL_" + block_str,
                                        mng_ptr_->thread_num_,
                                        mng_ptr_->output_file_ext_ + "_INTERNAL_" + block_str,
                                        mng_ptr_->inte_num_);
  internal_combine.process();

}

} /* namespace prot */
