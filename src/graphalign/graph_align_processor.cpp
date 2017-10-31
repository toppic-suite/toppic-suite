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
#include <string>
#include <vector>

#include "base/file_util.hpp"
#include "base/mod_util.hpp"
#include "base/threadpool.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "graph/proteo_graph_reader.hpp"
#include "graph/spec_graph_reader.hpp"
#include "graphalign/graph_align.hpp"
#include "graphalign/graph_align_processor.hpp"

namespace prot {

typedef std::shared_ptr<ThreadPool<PrsmXmlWriter> > PrsmXmlThreadPoolPtr;

std::function<void()> geneTask(FastaIndexReaderPtr reader_ptr,
                               GraphAlignMngPtr mng_ptr,
                               ModPtrVec var_mod_ptr_vec,
                               int spectrum_num, int idx) {
  return [reader_ptr, mng_ptr, var_mod_ptr_vec, spectrum_num, idx]() {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
    std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
    std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();

    std::string input_file_name
        = FileUtil::basename(sp_file_name) + "." + mng_ptr->input_file_ext_ + "_" + std::to_string(idx);
    SimplePrsmReader simple_prsm_reader(input_file_name);
    SimplePrsmStrPtr prsm_ptr = simple_prsm_reader.readOnePrsmStr();
    int group_spec_num = prsm_para_ptr->getGroupSpecNum();
    MsAlignReader sp_reader(sp_file_name, group_spec_num,
                            sp_para_ptr->getActivationPtr(),
                            sp_para_ptr->getSkipList());

    SpecGraphReader spec_reader(sp_file_name,
                                prsm_para_ptr->getGroupSpecNum(),
                                mng_ptr->convert_ratio_,
                                sp_para_ptr);

    PrsmXmlWriterPtr writer_ptr
        = std::make_shared<PrsmXmlWriter>(FileUtil::basename(sp_file_name) + "." + mng_ptr->output_file_ext_ + "_" + std::to_string(idx));
    SpectrumSetPtr spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0];
    ProteoAnnoPtr proteo_anno_ptr
        = std::make_shared<ProteoAnno>(prsm_para_ptr->getFixModPtrVec(),
                                       prsm_para_ptr->getProtModPtrVec(),
                                       var_mod_ptr_vec);

    int cnt = 0;
    while (spec_set_ptr != nullptr) {
      cnt += group_spec_num;
      if (spec_set_ptr->isValid()) {
        int spec_id = spec_set_ptr->getSpectrumId();
        std::vector<SimplePrsmStrPtr> selected_prsm_ptrs;
        while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
          selected_prsm_ptrs.push_back(prsm_ptr);
          prsm_ptr = simple_prsm_reader.readOnePrsmStr();
        }

        if (selected_prsm_ptrs.size() > 0) {
          SpecGraphPtrVec spec_ptr_vec
              = spec_reader.getNextSpecGraphPtrVec(spec_set_ptr, mng_ptr->prec_error_);
          for (size_t i = 0; i < selected_prsm_ptrs.size(); i++) {
            std::string seq_name = selected_prsm_ptrs[i]->getSeqName();
            std::string seq_desc = selected_prsm_ptrs[i]->getSeqDesc();
            std::vector<FastaSeqPtr> seq_ptr_vec = reader_ptr->readFastaSeqVec(seq_name, seq_desc);
            for (size_t j = 0; j < seq_ptr_vec.size(); j++) {
              for (size_t k = 0; k < spec_ptr_vec.size(); k++) {
                proteo_anno_ptr->anno(seq_ptr_vec[j]->getRawSeq(), seq_ptr_vec[j]->getSubSeqStart() == 0);
                MassGraphPtr graph_ptr = getMassGraphPtr(proteo_anno_ptr, mng_ptr->convert_ratio_);
                ProteoGraphPtr proteo_ptr = std::make_shared<ProteoGraph>(seq_ptr_vec[j],
                                                                          prsm_para_ptr->getFixModPtrVec(),
                                                                          graph_ptr,
                                                                          proteo_anno_ptr->isNme(),
                                                                          mng_ptr->convert_ratio_,
                                                                          mng_ptr->max_known_mods_,
                                                                          mng_ptr->getIntMaxPtmSumMass(),
                                                                          mng_ptr->proteo_graph_gap_,
                                                                          mng_ptr->var_ptm_in_gap_);
                GraphAlignPtr graph_align
                    = std::make_shared<GraphAlign>(mng_ptr, proteo_ptr, spec_ptr_vec[k]);
                graph_align->process();
                for (int shift = 0; shift <= mng_ptr->n_unknown_shift_; shift++) {
                  PrsmPtr prsm_ptr = graph_align->geneResult(shift);
                  if (prsm_ptr != nullptr) {
                    writer_ptr->write(prsm_ptr);
                  }
                }
                graph_align = nullptr;
              }
            }
          }
        }
      }
      if (idx == 0) {
        std::cout << std::flush << "Mass graph - processing " << cnt
            << " of " << spectrum_num << " spectra.\r";
      }
      spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0];
    }
    writer_ptr->close();
  };
}

void SimplePrsmFilter(SimplePrsmPtrVec & selected_prsm_ptrs) {
  if (selected_prsm_ptrs.size() == 0) return;
  if (selected_prsm_ptrs.size() == 1) {
    if (selected_prsm_ptrs[0]->getScore() < 10) selected_prsm_ptrs.clear();
    return;
  }
  std::sort(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(), SimplePrsm::cmpScoreDec);
  double scr = selected_prsm_ptrs[0]->getScore() * 2 / 3;
  selected_prsm_ptrs.erase(std::remove_if(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(),
                                          [scr] (const SimplePrsmPtr & p) {return p->getScore() < scr;}),
                           selected_prsm_ptrs.end());
  selected_prsm_ptrs.erase(std::remove_if(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(),
                                          [] (const SimplePrsmPtr & p) {return p->getScore() < 10;}),
                           selected_prsm_ptrs.end());
  selected_prsm_ptrs.erase(std::unique(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(),
                                       [] (const SimplePrsmPtr & a, const SimplePrsmPtr & b) {
                                         return a->getSpectrumScan() == b->getSpectrumScan()
                                           && a->getSeqName() == b->getSeqName();
                                       }),
                           selected_prsm_ptrs.end());
}

void GraphAlignProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  sp_para_ptr->prec_error_ = 0;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  LOG_DEBUG("Search db file name " << db_file_name);
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string var_mod_file_name = mng_ptr_->var_mod_file_name_;
  LOG_DEBUG("start reading " << var_mod_file_name);
  ModPtrVec var_mod_ptr_vec = ModUtil::readModTxt(var_mod_file_name)[2];
  LOG_DEBUG("end reading " << var_mod_file_name);

  int spectrum_num = MsAlignUtil::getSpNum(prsm_para_ptr->getSpectrumFileName());
  std::string input_file_name
      = FileUtil::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_;
  std::vector<std::shared_ptr<SimplePrsmXmlWriter> > simple_prsm_writer_vec;

  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    simple_prsm_writer_vec.push_back(std::make_shared<SimplePrsmXmlWriter>(input_file_name + "_" + std::to_string(i)));
  }

  SimplePrsmReader simple_prsm_reader(input_file_name);
  SimplePrsmPtr prsm_ptr = simple_prsm_reader.readOnePrsm();
  int cnt = 0;
  int spec_id = -1;
  SimplePrsmPtrVec selected_prsm_ptrs;
  while (prsm_ptr != nullptr) {
    if (prsm_ptr->getSpectrumId() == spec_id) {
      selected_prsm_ptrs.push_back(prsm_ptr); 
    } else {
      SimplePrsmFilter(selected_prsm_ptrs);
      cnt = cnt % mng_ptr_->thread_num_;
      simple_prsm_writer_vec[cnt]->write(selected_prsm_ptrs); 
      selected_prsm_ptrs.clear();
      cnt++;
      spec_id = prsm_ptr->getSpectrumId();
      selected_prsm_ptrs.push_back(prsm_ptr);
    }
    prsm_ptr = simple_prsm_reader.readOnePrsm(); 
  }
  simple_prsm_reader.close();

  SimplePrsmFilter(selected_prsm_ptrs);
  cnt = cnt % mng_ptr_->thread_num_;
  simple_prsm_writer_vec[cnt]->write(selected_prsm_ptrs); 

  for (size_t i = 0; i < simple_prsm_writer_vec.size(); i++) {
    simple_prsm_writer_vec[i]->close();
  }

  std::string output_file_name
      = FileUtil::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;

  FastaIndexReaderPtr reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);

  std::vector<ThreadPtr> thread_vec;
  for (int i = 1; i < mng_ptr_->thread_num_; i++) {
    ThreadPtr thread_ptr = std::make_shared<boost::thread>(geneTask(reader_ptr, mng_ptr_, var_mod_ptr_vec, spectrum_num, i));
    thread_vec.push_back(thread_ptr);
  }

  std::cout << std::flush << "Mass graph - processing 1 of " << spectrum_num << " spectra.\r";

  std::function<void()> mg_align_task = geneTask(reader_ptr, mng_ptr_, var_mod_ptr_vec, spectrum_num, 0);
  mg_align_task();

  for (size_t i = 0; i < thread_vec.size(); i++) {
    thread_vec[i]->join();
  }

  std::cout << std::flush << "Mass graph - processing " << spectrum_num
      << " of " << spectrum_num << " spectra." << std::endl;

  // combine result files
  std::vector<std::string> input_exts;
  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    std::string fname = mng_ptr_->output_file_ext_ + "_" + std::to_string(i);
    input_exts.push_back(fname);
  }

  int top_num = (mng_ptr_->n_unknown_shift_ + 1) * 4;
  PrsmStrCombinePtr combine_ptr
      = std::make_shared<PrsmStrCombine>(sp_file_name, input_exts,
                                         mng_ptr_->output_file_ext_, top_num);
  bool normalization = true;
  combine_ptr->process(normalization);
  combine_ptr = nullptr;
}
}  // namespace prot

