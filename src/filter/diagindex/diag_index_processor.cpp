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

#include <iostream>

#include "common/util/file_util.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/base_data.hpp"
#include "common/thread/simple_thread_pool.hpp"

#include "seq/proteoform.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform_factory.hpp"

#include "ms/spec/msalign_util.hpp"

#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_xml_writer_util.hpp"
#include "prsm/simple_prsm_str_merge.hpp"

#include "filter/diagindex/diag_index_processor.hpp"
#include "filter/diagindex/diag_index_file.hpp"

namespace toppic{

void DiagIndexProcessor::process() {
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::cout << "Generating Multi PTM index files --- started" << std::endl;

  std::vector<double> mod_mass_list;
  if (mng_ptr_->residueModFileName_ != "") {
    mod_mass_list = mod_util::getModMassVec(mod_util::readModTxt(mng_ptr_->residueModFileName_)[2]);
  }

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  int block_num = db_block_ptr_vec.size();
  int current_num = 1; //show how many files have been processed. n in the message "n of 5 files processed"..

  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    
    createIndexFiles(db_block_ptr_vec[i], db_block_ptr_vec.size(), mod_mass_list, block_num, &current_num);
  }
  pool_ptr->ShutDown();
  std::cout << "Generating Multi PTM index files --- finished" << std::endl;
}

void DiagIndexProcessor::createIndexFiles(DbBlockPtr block_ptr, int total_block_num,
                                       const std::vector<double> & mod_mass_list, int block_num, int *current_num) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;

  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
      + "_" + str_util::toString(block_ptr->getBlockIdx());

  ProteoformPtrVec raw_forms
      = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                        prsm_para_ptr->getFixModPtrVec());

  std::string block_str = str_util::toString(block_ptr->getBlockIdx());

  DiagIndexFilePtr index_file_ptr = std::make_shared<DiagIndexFile>(raw_forms, mng_ptr_, block_str);
  
  mng_ptr_->mutex_.lock();

    std::cout << "Multi PTM index files - processing " << *current_num << " of " << block_num << " files." << std::endl;
    *current_num = *current_num + 1;

    mng_ptr_->mutex_.unlock();

}


}