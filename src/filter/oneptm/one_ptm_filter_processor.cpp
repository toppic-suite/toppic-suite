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

#include <iomanip>
#include <iostream>

#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "common/base/mod_util.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

#include "spec/msalign_util.hpp"
#include "spec/spectrum_set.hpp"

#include "prsm/simple_prsm_xml_writer_set.hpp"
#include "prsm/simple_prsm_str_combine.hpp"

#include "filter/oneptm/one_ptm_filter_processor.hpp"
#include "filter/oneptm/mass_one_ptm_filter.hpp"

namespace toppic {


inline void filterBlock(const ProteoformPtrVec & raw_forms,
                        const std::string & block_str,
                        OnePtmFilterMngPtr mng_ptr,
                        const std::vector<double> & mod_mass_list) {
  MassOnePtmFilterPtr filter_ptr = std::make_shared<MassOnePtmFilter>(raw_forms, mng_ptr);
  int group_spec_num = mng_ptr->prsm_para_ptr_->getGroupSpecNum();
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(), 
                       group_spec_num,
                       prsm_para_ptr->getSpParaPtr()->getActivationPtr(),
                       prsm_para_ptr->getSpParaPtr()->getSkipList());
  std::string output_file_name = file_util::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr->output_file_ext_ + "_" + block_str;
  SimplePrsmXmlWriterSet writers(output_file_name);

  SpectrumSetPtrVec spec_set_vec = reader.getNextSpectrumSet(sp_para_ptr);

  while (spec_set_vec[0] != nullptr) {
    if (spec_set_vec[0]->isValid()) {
      if (mng_ptr->var_num_ == 0) {
        PrmMsPtrVec prm_ms_ptr_vec = spec_set_vec[0]->getMsTwoPtrVec();
        PrmMsPtrVec srm_ms_ptr_vec = spec_set_vec[0]->getSuffixMsTwoPtrVec();
        filter_ptr->computeBestMatch(prm_ms_ptr_vec, srm_ms_ptr_vec);
        writers.getCompleteWriterPtr()->write(filter_ptr->getCompMatchPtrs());
        writers.getPrefixWriterPtr()->write(filter_ptr->getPrefMatchPtrs());
        writers.getSuffixWriterPtr()->write(filter_ptr->getSuffMatchPtrs());
        writers.getInternalWriterPtr()->write(filter_ptr->getInternalMatchPtrs());
      } else {
        std::vector<double> mod_mass(3);
        for (size_t i = 0; i < mod_mass_list.size(); i++) {
          for (size_t k1 = 0; k1 < mod_mass.size(); k1++) {
            std::fill(mod_mass.begin(), mod_mass.end(), 0.0);
            mod_mass[k1] += mod_mass_list[i];
            PrmMsPtrVec prm_ms_ptr_vec = spec_set_vec[0]->getMsTwoPtrVec(sp_para_ptr, mod_mass);
            PrmMsPtrVec srm_ms_ptr_vec = spec_set_vec[0]->getSuffixMsTwoPtrVec(sp_para_ptr, mod_mass);
            filter_ptr->computeBestMatch(prm_ms_ptr_vec, srm_ms_ptr_vec);
            writers.getCompleteWriterPtr()->write(filter_ptr->getCompMatchPtrs());
            writers.getPrefixWriterPtr()->write(filter_ptr->getPrefMatchPtrs());
            writers.getSuffixWriterPtr()->write(filter_ptr->getSuffMatchPtrs());
            writers.getInternalWriterPtr()->write(filter_ptr->getInternalMatchPtrs());
          }
        }
      }
    }
    mng_ptr->cnt_++;
    double perc = mng_ptr->cnt_ * 100.0 / mng_ptr->n_spec_block_;
    std::stringstream msg;
    msg << std::flush << "One PTM filtering - processing " << std::setprecision(3) <<  perc << "%.     \r";
    std::cout << msg.str();
    spec_set_vec = reader.getNextSpectrumSet(sp_para_ptr);
  }
  reader.close();
  writers.close();
}

std::function<void()> geneTask(DbBlockPtr block_ptr,
                               const std::vector<double> &mod_mass_list, 
                               OnePtmFilterMngPtr mng_ptr) {
  return[block_ptr, mod_mass_list, mng_ptr] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
        + "_" + str_util::toString(block_ptr->getBlockIdx());
    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());
    std::string block_str = str_util::toString(block_ptr->getBlockIdx());
    filterBlock(raw_forms, block_str, mng_ptr, mod_mass_list);
  };
}

void OnePtmFilterProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::vector<double> mod_mass_list;
  if (mng_ptr_->residueModFileName_ != "") {
    mod_mass_list = mod_util::getModMassVec(mod_util::readModTxt(mng_ptr_->residueModFileName_)[2]);
  }

  // cnt_ is thread_safe 
  mng_ptr_->cnt_ = 0;
  int spec_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
  // n_spec_block = spec_num * block_num
  mng_ptr_->n_spec_block_ = spec_num * db_block_ptr_vec.size();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  int block_num = db_block_ptr_vec.size();
  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i], mod_mass_list, mng_ptr_));
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;

  std::cout << "One PTM filtering - combining blocks started." << std::endl;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_pref = mng_ptr_->output_file_ext_;
  SimplePrsmStrCombine::combineBlockResults(sp_file_name, input_pref, block_num,  
                                            mng_ptr_->comp_num_, mng_ptr_->pref_suff_num_, mng_ptr_->inte_num_ );
  std::cout << "One PTM filtering - combining blocks finished." << std::endl;
}

} /* namespace toppic */
