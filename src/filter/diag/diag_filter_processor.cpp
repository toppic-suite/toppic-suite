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
#include "common/thread/simple_thread_pool.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/factory/prm_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_xml_writer_util.hpp"
#include "prsm/simple_prsm_str_merge.hpp"

#include "filter/diag/diag_filter.hpp"
#include "filter/diag/diag_filter_processor.hpp"

namespace toppic {

inline void filterBlock(const ProteoformPtrVec & raw_forms,
                        int block_idx, 
                        DiagFilterMngPtr mng_ptr,
                        const std::vector<double> & mod_mass_list) {
  std::string block_str = str_util::toString(block_idx);

  DiagFilterPtr filter_ptr = std::make_shared<DiagFilter>(raw_forms, mng_ptr, block_str);

  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  int group_spec_num = mng_ptr->prsm_para_ptr_->getGroupSpecNum();

  SimpleMsAlignReaderPtr reader_ptr = std::make_shared<SimpleMsAlignReader>(sp_file_name,
                                                                            group_spec_num,
                                                                            sp_para_ptr->getActivationPtr());

  // init writer 
  std::string output_file_name = file_util::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr->output_file_ext_+"_"+ block_str; 

  SimplePrsmXmlWriterPtr writer_ptr = std::make_shared<SimplePrsmXmlWriter>(output_file_name);

  SpectrumSetPtr spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(reader_ptr, sp_para_ptr);

  while (spec_set_ptr != nullptr) {
    if (spec_set_ptr->isValid()) {
      if (mng_ptr->var_num_ == 0) {
        PrmMsPtrVec ms_ptr_vec = spec_set_ptr->getMsTwoPtrVec();
        SimplePrsmPtrVec match_ptrs = filter_ptr->getBestMatch(ms_ptr_vec);
        writer_ptr->write(match_ptrs);
      } 
      else {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        double prec_mono_mass = spec_set_ptr->getPrecMonoMass();
        std::vector<double> mod_mass(3);
        for (size_t i = 0; i < mod_mass_list.size(); i++) {
          for (size_t k1 = 0; k1 < mod_mass.size(); k1++) {
            std::fill(mod_mass.begin(), mod_mass.end(), 0.0);
            mod_mass[k1] += mod_mass_list[i];
            PrmMsPtrVec ms_ptr_vec = prm_ms_factory::geneMsTwoPtrVec(deconv_ms_ptr_vec,
                                                                     sp_para_ptr,
                                                                     prec_mono_mass, mod_mass);
            SimplePrsmPtrVec match_ptrs = filter_ptr->getBestMatch(ms_ptr_vec);
            writer_ptr->write(match_ptrs);
          }
        }
      }
    }
    mng_ptr->cnts_[block_idx] = mng_ptr->cnts_[block_idx] + group_spec_num;
    int cnt_sum = 0; 
    for (size_t i = 0; i < mng_ptr->cnts_.size(); i++) {
      cnt_sum = cnt_sum + mng_ptr->cnts_[i];
    }
    double perc = cnt_sum * 100.0 / mng_ptr->n_spec_block_;
    std::stringstream msg;
    msg << std::flush << "Multiple PTM filtering - processing " << std::setprecision(3) <<  perc << "%.     \r";
    mng_ptr->mutex_.lock();
    std::cout << msg.str();
    mng_ptr->mutex_.unlock();

    spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(reader_ptr, sp_para_ptr);
  }
  writer_ptr->close();
}

std::function<void()> geneTask(int block_idx, 
                               const std::vector<double> &mod_mass_list, 
                               DiagFilterMngPtr mng_ptr) {
  return[block_idx, mod_mass_list, mng_ptr] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
    std::string db_block_file_name = prsm_para_ptr->getOriDbName() + "_idx" 
      + file_util::getFileSeparator() + prsm_para_ptr->getSearchDbFileName()
      + "_" + str_util::toString(block_idx);
    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());
    filterBlock(raw_forms, block_idx, mng_ptr, mod_mass_list);
  };
}

void DiagFilterProcessor::process() {
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx" 
    + file_util::getFileSeparator() + mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::vector<double> mod_mass_list;
  if (mng_ptr_->residue_mod_file_name_ != "") {
    mod_mass_list 
        = mod_util::getModMassVec(mod_util::readModTxt(mng_ptr_->residue_mod_file_name_)[2]);
  }

  int spec_num = msalign_util::getSpNum(mng_ptr_->prsm_para_ptr_->getSpectrumFileName());
  mng_ptr_->n_spec_block_ = spec_num * db_block_ptr_vec.size();

  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  int block_num = db_block_ptr_vec.size();
  mng_ptr_->cnts_.resize(block_num, 0);

  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() == 0) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(10));
    }
    pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), mod_mass_list, mng_ptr_));
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;

  std::cout << "Multiple PTM filtering - combining blocks started." << std::endl;
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  SimplePrsmStrMergePtr merge_ptr
      = std::make_shared<SimplePrsmStrMerge>(sp_file_name, mng_ptr_->output_file_ext_,
                                             block_num, mng_ptr_->output_file_ext_,
                                             mng_ptr_->filter_result_num_);
  merge_ptr->process();
  merge_ptr = nullptr;
  //Remove temporary files
  file_util::cleanTempFiles(sp_file_name, mng_ptr_->output_file_ext_ + "_");
  std::cout << "Multiple PTM filtering - combining blocks finished." << std::endl;
}


}  // namespace toppic
