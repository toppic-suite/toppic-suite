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

#include <iomanip>
#include <iostream>

#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "common/base/mod_util.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

#include "ms/spec/msalign_util.hpp"
#include "ms/spec/spectrum_set.hpp"

#include "prsm/simple_prsm_xml_writer_set.hpp"
#include "prsm/simple_prsm_str_merge.hpp"

#include "filter/oneptm/one_ptm_filter_processor.hpp"
#include "filter/oneptm/mass_one_ptm_filter.hpp"
#include "filter/oneptm/mass_one_ptm_index_file.hpp"

#include "console/topindex_file_name.hpp"

namespace toppic {

inline void filterBlock(const ProteoformPtrVec & raw_forms,
                        int block_idx, 
                        OnePtmFilterMngPtr mng_ptr,
                        const std::vector<double> & mod_mass_list) {
  std::string block_str = str_util::toString(block_idx);
  MassOnePtmFilterPtr filter_ptr = std::make_shared<MassOnePtmFilter>(raw_forms, mng_ptr, block_str);
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

    mng_ptr->cnts_[block_idx] = mng_ptr->cnts_[block_idx] + 1;
    int cnt_sum = 0; 
    for (size_t i = 0; i < mng_ptr->cnts_.size(); i++) {
      cnt_sum = cnt_sum + mng_ptr->cnts_[i];
    }
    double perc = cnt_sum * 100.0 / mng_ptr->n_spec_block_;
    std::stringstream msg;
    msg << std::flush << "One PTM filtering - processing " << std::setprecision(3) <<  perc << "%.     \r";
    mng_ptr->mutex_.lock();
    std::cout << msg.str();
    mng_ptr->mutex_.unlock();
    spec_set_vec = reader.getNextSpectrumSet(sp_para_ptr);
  }
  reader.close();
  writers.close();
}

std::function<void()> geneTask(int block_idx, 
                               const std::vector<double> &mod_mass_list, 
                               OnePtmFilterMngPtr mng_ptr) {
  return[block_idx, mod_mass_list, mng_ptr] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
        + "_" + str_util::toString(block_idx);
    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());
    filterBlock(raw_forms, block_idx, mng_ptr, mod_mass_list);
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
  int spec_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
  // n_spec_block = spec_num * block_num
  mng_ptr_->n_spec_block_ = spec_num * db_block_ptr_vec.size();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  int block_num = db_block_ptr_vec.size();
  mng_ptr_->cnts_.resize(block_num, 0);
  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), mod_mass_list, mng_ptr_));
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;

  std::cout << "One PTM filtering - combining blocks started." << std::endl;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_pref = mng_ptr_->output_file_ext_;
  SimplePrsmStrMerge::mergeBlockResults(sp_file_name, input_pref, block_num,  
                                        mng_ptr_->comp_num_, mng_ptr_->pref_suff_num_, mng_ptr_->inte_num_ );
  std::cout << "One PTM filtering - combining blocks finished." << std::endl;
}
//below functions are used for generating index files

inline void createIndexFiles(const ProteoformPtrVec & raw_forms,
                        int block_idx, 
                        OnePtmFilterMngPtr mng_ptr,
                        const std::vector<double> & mod_mass_list, int block_num, int *current_num) {

    
    std::string block_str = str_util::toString(block_idx);
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;

    TopIndexFileName TopIndexFile;
    std::string parameters = TopIndexFile.gene_file_name(prsm_para_ptr);

    std::vector<std::string> file_vec{TopIndexFile.one_ptm_file_vec[0] + parameters + block_str, 
    TopIndexFile.one_ptm_file_vec[1] + parameters + block_str, 
    TopIndexFile.one_ptm_file_vec[2] + parameters + block_str, TopIndexFile.one_ptm_file_vec[3] + parameters + block_str};
    
    MassOnePtmIndexPtr filter_ptr = std::make_shared<MassOnePtmIndex>(raw_forms, mng_ptr, file_vec);
    
    mng_ptr->mutex_.lock();

    std::cout << "One PTM index files - processing " << *current_num << " of " << block_num << " files." << std::endl;
    *current_num = *current_num + 1;

    mng_ptr->mutex_.unlock();
}

std::function<void()> geneIndexTask(int block_idx, 
                               const std::vector<double> &mod_mass_list, 
                               OnePtmFilterMngPtr mng_ptr, int block_num, int *current_num) {
  return[block_idx, mod_mass_list, mng_ptr, block_num, current_num] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;

    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
        + "_" + str_util::toString(block_idx);
    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());
    createIndexFiles(raw_forms, block_idx, mng_ptr, mod_mass_list, block_num, current_num);
  };
}
void OnePtmFilterProcessor::index_process(){
  //for generating index files
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::cout << "Generating One PTM index files --- started" << std::endl;

  std::vector<double> mod_mass_list;
  if (mng_ptr_->residueModFileName_ != "") {
    mod_mass_list = mod_util::getModMassVec(mod_util::readModTxt(mng_ptr_->residueModFileName_)[2]);
  }

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  int block_num = db_block_ptr_vec.size();
  int current_num = 1; //show how many files have been processed. n in the message "n of 5 files processed"..

  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(geneIndexTask(db_block_ptr_vec[i]->getBlockIdx(), mod_mass_list, mng_ptr_, block_num, &current_num));
  }
  pool_ptr->ShutDown();
  std::cout << "Generating One PTM index files --- finished" << std::endl;
  //std::cout << std::endl;
}
} /* namespace toppic */
