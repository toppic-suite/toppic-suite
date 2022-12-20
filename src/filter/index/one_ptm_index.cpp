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

#include <boost/thread/mutex.hpp>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"

#include "seq/db_block.hpp"
#include "seq/proteoform_util.hpp"
#include "seq/proteoform_factory.hpp"

#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/index/one_ptm_index.hpp"

namespace toppic{

namespace one_ptm_index {

// serialization mutex.
boost::mutex serial_mutex;

void geneOnePtmIndexFile(const ProteoformPtrVec &proteo_ptrs,
                         OnePtmFilterMngPtr mng_ptr, std::vector<std::string> file_vec) {

  std::string dir_name = mng_ptr->prsm_para_ptr_->getOriDbName() + "_idx";
  std::vector<std::vector<double>> shift_2d
      = proteoform_util::getNTermShift2D(proteo_ptrs, 
                                         mng_ptr->prsm_para_ptr_->getProtModPtrVec());
  std::vector<std::vector<double>> n_term_acet_2d
      = proteoform_util::getNTermAcet2D(proteo_ptrs, 
                                        mng_ptr->prsm_para_ptr_->getProtModPtrVec());
  MassMatchPtr term_index_ptr 
      = mass_match_factory::getPrmTermMassMatchPtr(proteo_ptrs, shift_2d,
                                                   mng_ptr->max_proteoform_mass_,
                                                   mng_ptr->filter_scale_);
  {
    boost::unique_lock<boost::mutex> lock(serial_mutex);
    term_index_ptr->serializeMassMatch(file_vec[0], dir_name);
  }
  term_index_ptr = nullptr;

  MassMatchPtr diag_index_ptr 
      = mass_match_factory::getPrmDiagMassMatchPtr(proteo_ptrs,
                                                   mng_ptr->max_proteoform_mass_,
                                                   mng_ptr->filter_scale_);
  {
    boost::unique_lock<boost::mutex> lock(serial_mutex);
    diag_index_ptr->serializeMassMatch(file_vec[1], dir_name);
  }
  diag_index_ptr = nullptr;

  std::vector<std::vector<double>> rev_shift_2d;
  std::vector<double> shift_1d(1, 0);
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    rev_shift_2d.push_back(shift_1d);
  }
  MassMatchPtr rev_term_index_ptr 
      = mass_match_factory::getSrmTermMassMatchPtr(proteo_ptrs, 
                                                   rev_shift_2d, 
                                                   n_term_acet_2d,
                                                   mng_ptr->max_proteoform_mass_,
                                                   mng_ptr->filter_scale_);
  {
    boost::unique_lock<boost::mutex> lock(serial_mutex);
    rev_term_index_ptr->serializeMassMatch(file_vec[2], dir_name);
  }
  rev_term_index_ptr = nullptr;

  MassMatchPtr rev_diag_index_ptr 
      = mass_match_factory::getSrmDiagMassMatchPtr(proteo_ptrs, 
                                                   n_term_acet_2d,
                                                   mng_ptr->max_proteoform_mass_,
                                                   mng_ptr->filter_scale_);
  {
    boost::unique_lock<boost::mutex> lock(serial_mutex);
    rev_diag_index_ptr->serializeMassMatch(file_vec[3], dir_name); 
  }
  rev_diag_index_ptr = nullptr;
}

void createIndexFiles(const ProteoformPtrVec & raw_forms,
                      int block_idx, 
                      OnePtmFilterMngPtr mng_ptr) {
  std::string block_str = str_util::toString(block_idx);
  std::string parameters = mng_ptr->getIndexFilePara();
  std::vector<std::string> file_vec;
  for (size_t i = 0; i < mng_ptr->one_ptm_file_vec_.size(); i++){
    file_vec.push_back(mng_ptr->one_ptm_file_vec_[i] + parameters + block_str);
  }

  geneOnePtmIndexFile(raw_forms, mng_ptr, file_vec);
}

std::function<void()> geneIndexTask(int block_idx, 
                                    OnePtmFilterMngPtr mng_ptr) {
  return[block_idx, mng_ptr] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder()
      + "_" + str_util::toString(block_idx);
    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());
    createIndexFiles(raw_forms, block_idx, mng_ptr);
  };
}

void process(OnePtmFilterMngPtr mng_ptr) {
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::cout << "Generating one shift index files --- started" << std::endl;

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr->thread_num_);
  int block_num = db_block_ptr_vec.size();

  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() == 0) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    std::cout << "One shift index files - processing " << (i+1) 
        << " of " << block_num << " files." << std::endl;
    pool_ptr->Enqueue(geneIndexTask(db_block_ptr_vec[i]->getBlockIdx(), mng_ptr));
  }
  pool_ptr->ShutDown();
  std::cout << "Generating one shift index files --- finished" << std::endl;
}

}

}
