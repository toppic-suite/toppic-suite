// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include <utility>
#include <string>
#include <algorithm>
#include <vector>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/file_util.hpp"
#include "base/proteoform_factory.hpp"
#include "base/residue_util.hpp"
#include "base/fasta_reader.hpp"

#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "prsm/simple_prsm_table_writer.hpp"

#include "spec/msalign_util.hpp"

#include "suffix/db_file_handler.hpp"
#include "suffix/suffix_tree.hpp"

#include "tag_filter_processor.hpp"
#include "sequence_tag_graph.hpp"


namespace prot {

void TagFilterProcessor::process() {
  ModPtrVec fix_mod = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();
  ResiduePtrVec residue_list
      = residue_util::convertStrToResiduePtrVec("ARNDCEQGHIKMFPSTWYV", fix_mod);

  // generate mass list with variable ptms
  for (size_t i = 0; i < residue_list.size(); i++) {
    mass_list_.push_back(std::make_pair(residue_list[i]->getMass(),
                                        residue_list[i]->getAminoAcidPtr()->getOneLetter()));
  }

  std::sort(mass_list_.begin(), mass_list_.end(),
            [](const std::pair<double, std::string> & a,
               const std::pair<double, std::string> & b)->bool {
                 return a.first < b.first;
            });
  for (size_t i = 0; i < mass_list_.size(); i++) {
    LOG_DEBUG(mass_list_[i].first << " " << mass_list_[i].second)
  }
  processDB();
}

int TagFilterProcessor::compTagScore(const std::string & seq, const std::vector<int> & pos) {
  std::vector<int> hit_list(seq.length());
  std::fill(hit_list.begin(), hit_list.end(), 0);
  for (size_t i = 0; i < pos.size(); i++) {
    hit_list[pos[i]]++;
  }

  if (hit_list.size() < mng_ptr_->L) {
    return std::accumulate(hit_list.begin(), hit_list.end(), 0);
  }
  int max_scr = 0;
  size_t idx = 0;
  while (idx + mng_ptr_->L < hit_list.size()) {
    int tmp_scr = std::accumulate(hit_list.begin() + idx, hit_list.begin() + idx + mng_ptr_->L, 0);
    if (tmp_scr > max_scr) max_scr = tmp_scr;
    idx += 10;
  }
  return max_scr;
}

std::map<int, int> TagFilterProcessor::getHighScoreSeq(const std::vector<std::string> & tags,
                                                       suffix::SuffixTree *st,
                                                       suffix::ProteinDatabase *pd,
                                                       FastaIndexReaderPtr reader_ptr) {
  std::vector<suffix::SuffixPosition*> startPosList;
  for (size_t i = 0; i < tags.size(); i++) {
    LOG_DEBUG("searching tag " << tags[i]);
    std::vector<suffix::SuffixPosition*> tmp = st->search(tags[i]);
    startPosList.insert(startPosList.end(), tmp.begin(), tmp.end());
  }
  std::sort(startPosList.begin(), startPosList.end(),
            [] (suffix::SuffixPosition * a, suffix::SuffixPosition * b) {
              if (a->getSeqNum() == b->getSeqNum()) {
                return a->getPosInSeq() < b->getPosInSeq();
              }
              
              return a->getSeqNum() < b->getSeqNum();
            });

  std::map<int, int> seq_score_map;

  if (startPosList.empty()) return seq_score_map;

  int seq_num = startPosList[0]->getSeqNum();
  std::vector<int> pos;

  for (size_t i = 0; i < startPosList.size(); i++) {
    if (seq_num == startPosList[i]->getSeqNum()) {
      pos.push_back(startPosList[i]->getPosInSeq());
    } else {
      std::string seq
          = reader_ptr->readFastaSeq(pd->getProteinID(seq_num), pd->getProteinDesc(seq_num))->getRawSeq();
      seq_score_map[seq_num] = compTagScore(seq, pos);
      pos.clear();
      pos.push_back(startPosList[i]->getPosInSeq());
      seq_num = startPosList[i]->getSeqNum(); 
    }
  }
  return seq_score_map;
}

void TagFilterProcessor::processDB() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();

  LOG_DEBUG("db_file_name " << db_file_name);

  suffix::DatabaseFileHandler *df = new suffix::DatabaseFileHandler();
  suffix::ProteinDatabase *pd = df->loadDatabase(db_file_name);
  suffix::SuffixTree *st = new suffix::SuffixTree(pd->getSequence(), pd);

  FastaIndexReaderPtr reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);

  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr =  mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(),
                       group_spec_num,
                       sp_para_ptr->getActivationPtr(),
                       sp_para_ptr->getSkipList());

  std::string output_file_name
      = file_util::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;

  SimplePrsmXmlWriter writer(output_file_name);

  SpectrumSetPtr spec_set_ptr;

  while ((spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)[0]) != nullptr) {
    // basic spectrum graph
    DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
    LOG_DEBUG("Spec ID " <<  deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getId());
    if (spec_set_ptr->isValid()) {
      std::vector<PeakNodePtr> peaks;
      for (size_t k = 0; k < deconv_ms_ptr_vec[0]->size(); k++) {
        peaks.push_back(std::make_shared<PeakNode>(deconv_ms_ptr_vec[0]->getPeakPtr(k)->getMonoMass()));
      }

      std::sort(peaks.begin(), peaks.end(),
                [](PeakNodePtr a, PeakNodePtr b)->bool {
                return a->getMass() < b->getMass();
                });

      std::shared_ptr<SequenceTagGraph> seq_tag_graph = std::make_shared<SequenceTagGraph>();
      seq_tag_graph->setErrTolerance(mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo());
      seq_tag_graph->SetAminoAcidsArray(mass_list_);
      seq_tag_graph->SetPeakList(peaks);
      seq_tag_graph->SetNodeCount(peaks.size());
      seq_tag_graph->CollectSequenceTagGraphEdges();
      std::vector<std::string> seq_tag_set = seq_tag_graph->FindSequenceTags();

      LOG_DEBUG("seq_tag_set.size " << seq_tag_set.size());
      for (size_t i = 0; i < seq_tag_set.size(); i++) {
        LOG_DEBUG("Tag " << i << ": " << seq_tag_set[i]);
      }
      if (seq_tag_set.size() > 0) {
        std::map<int, int> seq_res = getHighScoreSeq(seq_tag_set, st, pd, reader_ptr);
        for (auto it = seq_res.begin(); it != seq_res.end(); it++) {
          SimplePrsmPtr simple_prsm
              = std::make_shared<SimplePrsm>(deconv_ms_ptr_vec[0]->getMsHeaderPtr(),
                                             group_spec_num,
                                             pd->getProteinID(it->first),
                                             pd->getProteinDesc(it->first),
                                             it->second);
          writer.write(simple_prsm);
        }
      }
    }
  }
  reader.close();
  writer.close();

  std::vector<std::string> input_exts;
  input_exts.push_back(mng_ptr_->output_file_ext_);

  SimplePrsmStrCombinePtr tag_combiner
      = std::make_shared<SimplePrsmStrCombine>(sp_file_name,
                                               input_exts,
                                               mng_ptr_->output_file_ext_ + std::to_string(mng_ptr_->top_num_),
                                               mng_ptr_->top_num_);

  tag_combiner->process();

  tag_combiner = nullptr;
}

}  // namespace prot
