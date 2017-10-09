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

#include <string>
#include <vector>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/db_block.hpp"
#include "base/proteoform_factory.hpp"
#include "base/file_util.hpp"
#include "base/mod_util.hpp"
#include "base/threadpool.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "diagfilter/mass_diag_filter.hpp"
#include "diagfilter/diag_filter_processor.hpp"

namespace prot {

typedef std::shared_ptr<ThreadPool<SimplePrsmXmlWriter>> SimplePrsmThreadPoolPtr;

std::function<void()> geneTask(MassDiagFilterPtr filter_ptr,
                               const PrmMsPtrVec & ms_ptr_vec,
                               SimplePrsmThreadPoolPtr  pool_ptr) {
  return [filter_ptr, ms_ptr_vec, pool_ptr]() {
    SimplePrsmPtrVec match_ptrs = filter_ptr->getBestMatch(ms_ptr_vec);
    boost::thread::id thread_id = boost::this_thread::get_id();
    std::shared_ptr<SimplePrsmXmlWriter> writer_ptr = pool_ptr->getWriter(thread_id);
    writer_ptr->write(match_ptrs);
  };
}

void DiagFilterProcessor::process() {
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::vector<double> mod_mass_list;
  if (mng_ptr_->residueModFileName_ != "") {
    mod_mass_list = ModUtil::getModMassVec(ModUtil::readModTxt(mng_ptr_->residueModFileName_)[2]);
  }

  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    std::cout << "Diagonal filtering - block " << (i + 1) << " out of "
        << db_block_ptr_vec.size() << " started." << std::endl;
    processBlock(db_block_ptr_vec[i], db_block_ptr_vec.size(), mod_mass_list);
    std::cout << "Diagonal filtering - block " << (i + 1) << " finished. " << std::endl;
  }

  std::cout << "Diagonal filtering - combining blocks started." << std::endl;

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int block_num = db_block_ptr_vec.size();

  SimplePrsmStrCombinePtr combine_ptr
      = std::make_shared<SimplePrsmStrCombine>(sp_file_name, mng_ptr_->output_file_ext_,
                                               block_num, mng_ptr_->output_file_ext_,
                                               mng_ptr_->filter_result_num_);
  combine_ptr->process();
  combine_ptr = nullptr;
  std::cout << "Diagonal filtering - combining blocks finished." << std::endl;
}

void DiagFilterProcessor::processBlock(DbBlockPtr block_ptr, int total_block_num,
                                       const std::vector<double> & mod_mass_list) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
      + "_" + std::to_string(block_ptr->getBlockIdx());
  ProteoformPtrVec raw_forms
      = ProteoformFactory::readFastaToProteoformPtrVec(db_block_file_name,
                                                       prsm_para_ptr->getFixModPtrVec());
  MassDiagFilterPtr filter_ptr = std::make_shared<MassDiagFilter>(raw_forms, mng_ptr_);

  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr =  mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(),
                       group_spec_num,
                       sp_para_ptr->getActivationPtr(),
                       sp_para_ptr->getSkipList());

  std::string output_file_name = FileUtil::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_+"_"+ std::to_string(block_ptr->getBlockIdx());

  SimplePrsmThreadPoolPtr pool_ptr
      = std::make_shared<ThreadPool<SimplePrsmXmlWriter>>(mng_ptr_->thread_num_, output_file_name);

  SpectrumSetPtr spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)[0];
  int spectrum_num = MsAlignUtil::getSpNum(prsm_para_ptr->getSpectrumFileName());
  int cnt = 0;
  while (spec_set_ptr != nullptr) {
    cnt += group_spec_num;
    if (spec_set_ptr->isValid()) {
      if (mng_ptr_->var_num_ == 0) {
        PrmMsPtrVec ms_ptr_vec = spec_set_ptr->getMsTwoPtrVec();
        while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
          boost::this_thread::sleep(boost::posix_time::milliseconds(100));
        }
        pool_ptr->Enqueue(geneTask(filter_ptr, ms_ptr_vec, pool_ptr));
      } else {
        for (size_t i = 0; i < mod_mass_list.size(); i++) {
          for (size_t k1 = 0; k1 < sp_para_ptr->mod_mass_.size(); k1++) {
            std::fill(sp_para_ptr->mod_mass_.begin(), sp_para_ptr->mod_mass_.end(), 0.0);
            sp_para_ptr->mod_mass_[k1] += mod_mass_list[i];
            PrmMsPtrVec ms_ptr_vec = spec_set_ptr->getMsTwoPtrVec(sp_para_ptr);
            while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
              boost::this_thread::sleep(boost::posix_time::milliseconds(100));
            }
            pool_ptr->Enqueue(geneTask(filter_ptr, ms_ptr_vec, pool_ptr));
          }
        }
      }
    }
    std::cout << std::flush << "Diagonal filtering - processing " << cnt
        << " of " << spectrum_num << " spectra.\r";
    spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)[0];
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;
  reader.close();
  std::vector<std::string> input_exts;
  std::string cur_output_ext =
      mng_ptr_->output_file_ext_ + "_" + std::to_string(block_ptr->getBlockIdx());
  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    std::string fname = cur_output_ext + "_" + std::to_string(i);
    input_exts.push_back(fname);
  }

  SimplePrsmStrCombinePtr combine_ptr
      = std::make_shared<SimplePrsmStrCombine>(mng_ptr_->prsm_para_ptr_->getSpectrumFileName(),
                                               input_exts, cur_output_ext, INT_MAX);
  combine_ptr->process();
  combine_ptr = nullptr;
}

}  // namespace prot
