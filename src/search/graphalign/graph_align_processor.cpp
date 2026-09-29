// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "search/graphalign/graph_align_processor.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "common/base/mod_util.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "common/util/mem_check.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_xml_writer_util.hpp"
#include "search/graph/proteo_graph_reader.hpp"
#include "search/graph/spec_graph_reader.hpp"
#include "search/graphalign/graph_align.hpp"
#include "seq/fasta_sub_util.hpp"

namespace toppic {

// The alignments of one candidate protein of a spectrum: every sub-sequence
// of the protein against every spectrum graph of the spectrum, and the PrSM
// of each unexpected shift number, as the XML strings of the writer of the
// candidate's file.
void alignCandidate(const GraphAlignMngPtr& mng_ptr,
                    const SimplePrsmStrPtr& cand_ptr,
                    const SpecGraphPtrVec& spec_ptr_vec,
                    const FastaIndexReaderPtr& reader_ptr,
                    const ProteoAnnoPtr& proteo_anno_ptr,
                    const PrsmXmlWriterPtr& writer_ptr,
                    std::vector<std::string>* results) {
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  std::string seq_name = cand_ptr->getSeqName();
  std::string seq_desc = cand_ptr->getSeqDesc();
  FastaSeqPtr seq_ptr = reader_ptr->readFastaSeq(seq_name, seq_desc);
  std::vector<FastaSubSeqPtr> seq_ptr_vec = fasta_sub_util::breakSeq(seq_ptr);

  for (size_t j = 0; j < seq_ptr_vec.size(); j++) {
    for (size_t k = 0; k < spec_ptr_vec.size(); k++) {
      proteo_anno_ptr->anno(seq_ptr_vec[j]->getRawSeq(),
                            seq_ptr_vec[j]->isNTerm());
      MassGraphPtr graph_ptr =
          getMassGraphPtr(proteo_anno_ptr, mng_ptr->convert_ratio_);
      ProteoGraphPtr proteo_ptr = std::make_shared<ProteoGraph>(
          seq_ptr_vec[j], prsm_para_ptr->getFixModPtrVec(), graph_ptr,
          proteo_anno_ptr->isNme(), mng_ptr->convert_ratio_,
          mng_ptr->max_known_mods_, mng_ptr->getIntMaxPtmSumMass(),
          mng_ptr->proteo_graph_gap_, mng_ptr->var_ptm_in_gap_,
          &spec_ptr_vec[k]->getDistCover());
      GraphAlignPtr graph_align = std::make_shared<GraphAlign>(
          mng_ptr, proteo_ptr, spec_ptr_vec[k], seq_ptr);
      graph_align->process();
      for (int shift = 0; shift <= mng_ptr->n_unknown_shift_; shift++) {
        PrsmPtr prsm_ptr = graph_align->geneResult(shift);
        if (prsm_ptr != nullptr) {
          results->push_back(writer_ptr->toXmlString(prsm_ptr));
        }
      }
    }
  }
}

void SimplePrsmFilter(SimplePrsmPtrVec& selected_prsm_ptrs) {
  if (selected_prsm_ptrs.size() == 0) return;
  if (selected_prsm_ptrs.size() == 1) {
    if (selected_prsm_ptrs[0]->getScore() < 10) selected_prsm_ptrs.clear();
    return;
  }
  std::sort(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(),
            SimplePrsm::cmpScoreDecSeqInc);
  selected_prsm_ptrs.erase(
      std::remove_if(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(),
                     [](const SimplePrsmPtr& p) { return p->getScore() < 10; }),
      selected_prsm_ptrs.end());
  selected_prsm_ptrs.erase(
      std::unique(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(),
                  [](const SimplePrsmPtr& a, const SimplePrsmPtr& b) {
                    return a->getSpectrumScan() == b->getSpectrumScan() &&
                           a->getSeqName() == b->getSeqName();
                  }),
      selected_prsm_ptrs.end());
}

void GraphAlignProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();

  LOG_DEBUG("Search db file name " << db_file_name);
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string var_mod_file_name = mng_ptr_->var_mod_file_name_;
  LOG_DEBUG("start reading " << var_mod_file_name);
  ModPtrVec var_mod_ptr_vec = mod_util::readModTxt(var_mod_file_name)[2];
  LOG_DEBUG("end reading " << var_mod_file_name);

  int spectrum_num =
      msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
  std::string input_file_name =
      file_util::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_;

  SimplePrsmReaderPtr simple_prsm_reader =
      std::make_shared<toppic::SimplePrsmReader>(input_file_name);
  std::shared_ptr<SimplePrsmXmlWriter> graph_filter_writer =
      std::make_shared<SimplePrsmXmlWriter>(file_util::basename(sp_file_name) +
                                            ".topmg_graph");
  SimplePrsmPtr prsm_ptr = simple_prsm_reader->readOnePrsm();
  int spec_id = -1;
  SimplePrsmPtrVec selected_prsm_ptrs;
  // only the spectra of the list are searched: the candidates of the others
  // are dropped here, so neither the alignment nor anything after it sees
  // them
  std::set<int> searched_spec_ids;
  std::set<int> skipped_spec_ids;
  while (prsm_ptr != nullptr) {
    if (!mng_ptr_->search_all_ &&
        mng_ptr_->search_spec_ids_.count(prsm_ptr->getSpectrumId()) == 0) {
      skipped_spec_ids.insert(prsm_ptr->getSpectrumId());
      prsm_ptr = simple_prsm_reader->readOnePrsm();
      continue;
    }
    searched_spec_ids.insert(prsm_ptr->getSpectrumId());
    if (prsm_ptr->getSpectrumId() == spec_id) {
      selected_prsm_ptrs.push_back(prsm_ptr);
    } else {
      SimplePrsmFilter(selected_prsm_ptrs);
      graph_filter_writer->write(selected_prsm_ptrs);
      selected_prsm_ptrs.clear();
      spec_id = prsm_ptr->getSpectrumId();
      selected_prsm_ptrs.push_back(prsm_ptr);
    }
    prsm_ptr = simple_prsm_reader->readOnePrsm();
  }
  simple_prsm_reader->close();
  SimplePrsmFilter(selected_prsm_ptrs);
  graph_filter_writer->write(selected_prsm_ptrs);
  graph_filter_writer->close();
  if (!mng_ptr_->search_all_) {
    std::cout << "Graph alignment - " << searched_spec_ids.size()
              << " of the spectra to search have candidates, "
              << skipped_spec_ids.size() << " other spectra skipped."
              << std::endl;
  }

  // The candidates of the .topmg_graph file are aligned by whichever
  // thread is free. Candidate c used to be aligned by thread c % thread_num,
  // which took the candidates of its own file in turn and wrote their PrSMs
  // to its own output file, and the merge below reads the PrSMs of a
  // spectrum file by file; the PrSMs are therefore kept per candidate and
  // written to the same files in the same order once every candidate is
  // aligned, which leaves the merged file as it was. A thread's turn could
  // not get past a candidate whose spectrum it never met, so such a
  // candidate stops its file here too.
  //
  // The files are still thread_num_ of them, but the thread pool may be
  // smaller: every thread holds an alignment's DP table and consistent
  // pairs, and mem_check lowers the number of threads to the memory
  // available (toppic_test experiment 77). The output does not depend on
  // the pool's size, only on the number of files.
  int thread_num = mng_ptr_->thread_num_;
  int layer_num = (mng_ptr_->n_unknown_shift_ + 1) *
                  (mng_ptr_->max_known_mods_ + 1);
  int pool_thread_num =
      mem_check::getGraphAlignThreadNum(thread_num, layer_num);
  SimplePrsmStrPtrVec cands;
  {
    SimplePrsmReader cand_reader(file_util::basename(sp_file_name) +
                                 ".topmg_graph");
    SimplePrsmStrPtr cand_ptr = cand_reader.readOnePrsmStr();
    while (cand_ptr != nullptr) {
      cands.push_back(cand_ptr);
      cand_ptr = cand_reader.readOnePrsmStr();
    }
    cand_reader.close();
  }

  PrsmXmlWriterPtrVec writer_ptrs;
  std::vector<FastaIndexReaderPtr> reader_ptrs;
  std::vector<ProteoAnnoPtr> proteo_anno_ptrs;
  for (int i = 0; i < thread_num; i++) {
    writer_ptrs.push_back(std::make_shared<PrsmXmlWriter>(
        file_util::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_ +
        "_" + std::to_string(i)));
  }
  for (int i = 0; i < pool_thread_num; i++) {
    // what a thread changes or holds open: one per thread of the pool
    reader_ptrs.push_back(std::make_shared<FastaIndexReader>(db_file_name));
    proteo_anno_ptrs.push_back(std::make_shared<ProteoAnno>(
        prsm_para_ptr->getFixModPtrVec(), prsm_para_ptr->getProtModPtrVec(),
        var_mod_ptr_vec));
  }
  std::vector<std::vector<std::string>> results(cands.size());
  std::vector<bool> file_stopped(thread_num, false);

  SimpleThreadPoolPtr pool_ptr =
      std::make_shared<SimpleThreadPool>(pool_thread_num);
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReaderPtr ms_reader_ptr = std::make_shared<MsAlignReader>(
      sp_file_name, group_spec_num, sp_para_ptr->getActivationPtr(),
      sp_para_ptr->getEnvCnnCutoff(), sp_para_ptr->getFragDedupPpo());
  SpecGraphReader spec_reader(sp_file_name, group_spec_num,
                              mng_ptr_->convert_ratio_, sp_para_ptr);
  size_t next_cand = 0;
  int cnt = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  while (deconv_ms_ptr_vec.size() > 0) {
    if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->containsPrec()) {
      double prec_mono_mass =
          deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass();
      SpectrumSetPtr spec_set_ptr = spectrum_set_factory::geneSpectrumSetPtr(
          deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
      if (spec_set_ptr->isValid()) {
        int spec_id = spec_set_ptr->getSpectrumId();
        while (next_cand < cands.size() &&
               cands[next_cand]->getSpectrumId() < spec_id) {
          file_stopped[next_cand % thread_num] = true;
          next_cand++;
        }
        size_t end_cand = next_cand;
        while (end_cand < cands.size() &&
               cands[end_cand]->getSpectrumId() == spec_id) {
          end_cand++;
        }
        if (end_cand > next_cand) {
          // the spectrum graphs are built once and shared by its candidates
          std::shared_ptr<SpecGraphPtrVec> spec_ptr_vec =
              std::make_shared<SpecGraphPtrVec>(
                  spec_reader.getNextSpecGraphPtrVec(spec_set_ptr,
                                                     mng_ptr_->prec_error_));
          // the proteoform distances the alignments can pair with each
          for (size_t k = 0; k < spec_ptr_vec->size(); k++) {
            (*spec_ptr_vec)[k]->compDistCover(
                mng_ptr_->getIntMinConsistentDist(),
                mng_ptr_->getIntTolerance());
          }
          for (size_t c = next_cand; c < end_cand; c++) {
            if (file_stopped[c % thread_num]) {
              continue;
            }
            // a short queue, so that only a few spectra are held
            while (pool_ptr->getQueueSize() >=
                   static_cast<size_t>(2 * pool_thread_num)) {
              std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            GraphAlignMngPtr mng_ptr = mng_ptr_;
            SimplePrsmStrPtr cand_ptr = cands[c];
            PrsmXmlWriterPtr writer_ptr = writer_ptrs[c % thread_num];
            std::vector<std::string>* cand_results = &results[c];
            SimpleThreadPool* pool = pool_ptr.get();
            pool_ptr->enqueue([mng_ptr, cand_ptr, spec_ptr_vec, writer_ptr,
                               cand_results, pool, &reader_ptrs,
                               &proteo_anno_ptrs]() {
              int id = pool->getId(std::this_thread::get_id());
              alignCandidate(mng_ptr, cand_ptr, *spec_ptr_vec,
                             reader_ptrs[id], proteo_anno_ptrs[id],
                             writer_ptr, cand_results);
            });
          }
          next_cand = end_cand;
        }
      }
    }
    cnt += group_spec_num;
    std::cout << std::flush << "Mass graph alignment - processing " << cnt
              << " of " << spectrum_num << " spectra.\r";
    deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  }
  pool_ptr->shutDown();
  std::cout << std::flush << "Mass graph alignment - processing "
            << spectrum_num << " of " << spectrum_num << " spectra."
            << std::endl;

  for (size_t c = 0; c < results.size(); c++) {
    for (size_t k = 0; k < results[c].size(); k++) {
      writer_ptrs[c % thread_num]->writeXmlString(results[c][k]);
    }
  }
  for (int i = 0; i < thread_num; i++) {
    writer_ptrs[i]->close();
  }

  // combine result files
  std::vector<std::string> input_exts;
  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    std::string fname = mng_ptr_->output_file_ext_ + "_" + std::to_string(i);
    input_exts.push_back(fname);
  }

  int top_num = (mng_ptr_->n_unknown_shift_ + 1) * 4;
  bool norm = true;
  bool remove_dup = true;
  PrsmStrMergePtr merge_ptr = std::make_shared<PrsmStrMerge>(
      sp_file_name, input_exts, mng_ptr_->output_file_ext_, top_num, norm,
      remove_dup);
  merge_ptr->process();
  merge_ptr = nullptr;

  // remove temporary files
  for (int t = 0; t < mng_ptr_->thread_num_; t++) {
    file_util::cleanTempFiles(
        sp_file_name, mng_ptr_->output_file_ext_ + "_" + std::to_string(t));
  }
}
}  // namespace toppic
