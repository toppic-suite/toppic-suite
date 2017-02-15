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


#include <iostream>
#include <fstream>

#include "base/file_util.hpp"
#include "base/mod_util.hpp"
#include "base/threadpool.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "graph/graph_util.hpp"
#include "graph/proteo_graph_reader.hpp"
#include "graph/spec_graph_reader.hpp"
#include "graphalign/graph_align.hpp"
#include "graphalign/graph_align_processor.hpp"

namespace prot {

typedef std::shared_ptr<ThreadPool<PrsmXmlWriter>> PrsmXmlThreadPoolPtr;

std::function<void()> geneTask(GraphAlignMngPtr mng_ptr,
                               ProteoGraphPtr proteo_ptr,
                               SpecGraphPtr spec_ptr,
                               PrsmXmlThreadPoolPtr pool_ptr) {
  return [mng_ptr, proteo_ptr, spec_ptr, pool_ptr]() {
    GraphAlignPtr graph_align 
        = std::make_shared<GraphAlign>(mng_ptr, proteo_ptr, spec_ptr);
    graph_align->process();
    boost::thread::id thread_id = boost::this_thread::get_id();
    PrsmXmlWriterPtr writer_ptr = pool_ptr->getWriter(thread_id);
    for (int shift = 0; shift <= mng_ptr->n_unknown_shift_; shift++) {
      PrsmPtr prsm_ptr = graph_align->geneResult(shift);
      if (prsm_ptr != nullptr) {
        writer_ptr->write(prsm_ptr);
      }
    }
    graph_align = nullptr;
  };
}

void GraphAlignProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  LOG_DEBUG("Search db file name " << db_file_name);
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string var_mod_file_name = mng_ptr_->var_mod_file_name_;
  LOG_DEBUG("start reading " << var_mod_file_name);
  ModPtrVec var_mod_ptr_vec = ModUtil::readModTxt(var_mod_file_name)[2];
  LOG_DEBUG("end reading " << var_mod_file_name);

  ProteoGraphReader reader(db_file_name, 
                           prsm_para_ptr->getFixModPtrVec(),
                           prsm_para_ptr->getProtModPtrVec(),
                           var_mod_ptr_vec,
                           mng_ptr_->convert_ratio_,
                           mng_ptr_->max_known_mods_,
                           mng_ptr_->getIntMaxPtmSumMass(),
                           mng_ptr_->proteo_graph_gap_);
  LOG_DEBUG("init reader complete");

  int spectrum_num = MsAlignUtil::getSpNum (prsm_para_ptr->getSpectrumFileName());
  std::string input_file_name = FileUtil::basename(sp_file_name)+ "." + mng_ptr_->input_file_ext_;
  SimplePrsmReader simple_prsm_reader(input_file_name);
  SimplePrsmPtr prsm_ptr = simple_prsm_reader.readOnePrsm();
  std::string output_file_name = FileUtil::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          sp_para_ptr->getActivationPtr());

  SpecGraphReader spec_reader(sp_file_name, 
                              prsm_para_ptr->getGroupSpecNum(),
                              mng_ptr_->convert_ratio_,
                              sp_para_ptr);

  PrsmXmlThreadPoolPtr pool_ptr =
      std::make_shared<ThreadPool<PrsmXmlWriter>>(mng_ptr_->thread_num_, output_file_name);

  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;

  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    cnt += group_spec_num;
    if(spec_set_ptr->isValid()){
      int spec_id = spec_set_ptr->getSpecId();
      SimplePrsmPtrVec selected_prsm_ptrs;
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        if (prsm_ptr->getScore() > 4)
          selected_prsm_ptrs.push_back(prsm_ptr);
        prsm_ptr = simple_prsm_reader.readOnePrsm();
      }

      if (selected_prsm_ptrs.size() > 0) {
        SimplePrsmPtrVec simple_prsm_ptrs = 
            SimplePrsmUtil::getUniqueMatches(selected_prsm_ptrs);
        FastaIndexReaderPtr reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);
        for (size_t i = 0; i < simple_prsm_ptrs.size(); i++) {
          std::string seq_name = simple_prsm_ptrs[i]->getSeqName();
          std::string seq_desc = simple_prsm_ptrs[i]->getSeqDesc();
          FastaSeqPtr seq_ptr = reader_ptr->readFastaSeq(seq_name, seq_desc);
          ProteoGraphPtr proteo_ptr = reader.getProteoGraphPtrBySeq(seq_ptr);
          SpecGraphPtrVec spec_ptr_vec = spec_reader.getNextSpecGraphPtrVec(spec_set_ptr, mng_ptr_->prec_error_);
          for (size_t k = 0; k < spec_ptr_vec.size(); k++) {
            while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
              boost::this_thread::sleep(boost::posix_time::milliseconds(100));
            }
            pool_ptr->Enqueue(geneTask(mng_ptr_, proteo_ptr,  spec_ptr_vec[k], pool_ptr));
          }
        }
      }
    }
    std::cout << std::flush <<  "Mass graph - processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";
  }

  pool_ptr->ShutDown();
  std::cout << std::endl;

  // combine result files
  std::vector<std::string> input_exts;
  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    std::string fname = mng_ptr_->output_file_ext_ + "_" + std::to_string(i); 
    input_exts.push_back(fname);
  }
  int top_num = 1;
  PrsmStrCombinePtr combine_ptr
      = std::make_shared<PrsmStrCombine>(sp_file_name, input_exts,
                                         mng_ptr_->output_file_ext_, top_num);
  bool normalization = true;
  combine_ptr->process(normalization);
  combine_ptr = nullptr;

}

}

