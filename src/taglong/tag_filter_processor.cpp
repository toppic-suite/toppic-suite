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


#include <regex>
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

#include "tag_filter_processor.hpp"
#include "peak_node.hpp"


namespace prot {

void TagFilterProcessor::process() {
  ModPtrVec fix_mod = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();
  ResiduePtrVec residue_list
      = residue_util::convertStrToResiduePtrVec("ARNDCEQGHILKMFPSTWYVUO", fix_mod);

  // generate mass list with variable ptms
  for (size_t i = 0; i < residue_list.size(); i++) {
    mass_list_.push_back(std::make_pair(residue_list[i]->getMass(),
                                        residue_list[i]->getAcidPtr()->getOneLetter()));
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

std::vector<double> TagFilterProcessor::getEdgeLimits(PeakNodePtr peak, PeakNodePtr next) {
  double diff = next->getMass() - peak->getMass();
  std::vector<double> limits;
  double firstMass = peak->getMass();
  double secondMass = next->getMass();

  // larger error tolerance needed for combined spectrum graph
  double error =  (firstMass + secondMass) * mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo() / 2;
  limits.push_back(diff - error);
  limits.push_back(diff + error);
  return limits;
}

void TagFilterProcessor::generateGapEdges(std::vector<PeakNodePtr> & peaks) {
  for (size_t i = 0; i < peaks.size(); i++) {
    peaks[i]->clearEdges();
  }

  for (size_t i = 0; i < peaks.size() - 1; i++) {
    for (size_t j = i + 1; j < peaks.size(); j++) {
      PeakNodePtr next = peaks[j];
      std::vector<double> limits = getEdgeLimits(peaks[i], next);
      for (size_t k = 0; k < mass_list_.size(); k++) {
        if (limits[0] <= mass_list_[k].first && limits[1] >= mass_list_[k].first) {
          peaks[i]->addNext(next);
          break;
        }
      }
    }
  }
}

std::vector<std::vector<PeakNodePtr> > getComponentsFromGraph(std::vector<PeakNodePtr> peaks) {
  for (size_t i = 0; i < peaks.size(); i++) {
    peaks[i]->setComponenetId(i);
  }

  bool done = false;

  do {
    done = true;
    for (size_t i = 0; i < peaks.size(); i++) {
      if (peaks[i]->updateComponentId()) {
        done = false;
      }
    }
  } while (!done);

  std::vector<bool> componentDone(peaks.size());
  std::fill(componentDone.begin(), componentDone.end(), false);

  std::vector<std::vector<PeakNodePtr> > components;

  for (size_t i = 0 ; i < peaks.size(); i++) {
    PeakNodePtr p = peaks[i];
    int componentId = p->getComponenetId();
    if (!componentDone[componentId]) {
      std::vector<PeakNodePtr> component;
      for (size_t j = 0; j < peaks.size(); j++) {
        if (peaks[j]->getComponenetId() == componentId) {
          component.push_back(peaks[j]);
        }
      }

      if (component.size() > 1) {
        components.push_back(component);
        componentDone[componentId] = true;
      }
    }
  }

  return components;
}

std::vector<PeakNodePtr> findBestTag(PeakNodePtr peak, std::vector<PeakNodePtr> best,
                                     int len, std::vector<PeakNodePtr> & prefix) {
  if (peak->getMaxPrefix() >= len) {
    return best;
  }

  prefix[len] = peak;

  if (len >= static_cast<int>(best.size())) {
    best.resize(len + 1);
    for (int i = 0 ; i <= len; i++) {
      best[i] = prefix[i];
    }
  }

  for (size_t i = 0 ; i < peak->getNext().size(); i++) {
    PeakNodePtr next = peak->getNext()[i];
    best = findBestTag(next, best, len + 1, prefix);
  }

  peak->setMaxPrefix(len);

  return best;
}

std::vector<PeakNodePtr> findBestTag(std::vector<PeakNodePtr> peaks) {
  std::vector<PeakNodePtr> bestTag;

  for (size_t i = 0; i < peaks.size(); i++) {
    peaks[i]->setMaxPrefix(-1);
  }

  std::vector<PeakNodePtr> prefix;
  for (int i = 0; i < 500; i++) {
    prefix.push_back(std::make_shared<PeakNode>(0.0));
  }

  for (size_t i = 0; i < peaks.size(); i++) {
    bestTag = findBestTag(peaks[i], bestTag, 0, prefix);
  }
  return bestTag;
}

std::vector<std::string> TagFilterProcessor::getTags(std::vector<std::vector<PeakNodePtr> > componentsFromGraph) {
  std::vector<std::vector<std::string> > tags;

  for (size_t i = 0; i < componentsFromGraph.size(); i++) {
    std::vector<PeakNodePtr> tagPeaks = findBestTag(componentsFromGraph[i]);
    std::vector<std::string> superTag;
    if (tagPeaks.size() < mng_ptr_->tag_min_len_ + 1) {
      continue;
    }
    for (size_t j = 0; j < tagPeaks.size() - 1; j++) {
      std::vector<double> limits = getEdgeLimits(tagPeaks[j], tagPeaks[j + 1]);
      for (size_t k = 0; k < mass_list_.size(); k++) {
        if (limits[0] <= mass_list_[k].first && limits[1] >= mass_list_[k].first) {
          superTag.push_back(mass_list_[k].second);
          break;
        }
      }
    }
    tags.push_back(superTag);
  }

  std::vector<std::string> break_tags;
  for (size_t i = 0; i < tags.size(); i++) {
    if (tags[i].size() < mng_ptr_->tag_min_len_) {
      std::abort();
    }

    if (tags[i].size() == mng_ptr_->tag_min_len_) {
      std::string tmp = "";
      std::string tmp2 = "";

      for (size_t k = 0; k < tags[i].size(); k++) {
        tmp = tmp + tags[i][k];
        tmp2 = tags[i][k] + tmp2;
      }
      break_tags.push_back(tmp);
      break_tags.push_back(tmp2);
    } else {
      for (size_t j = 0; j < tags[i].size() - mng_ptr_->tag_min_len_ + 1; j++) {
        std::string tmp = "";
        std::string tmp2 = "";
        for (size_t k = j; k < j + mng_ptr_->tag_min_len_; k++) {
          tmp = tmp + tags[i][k];
          tmp2 = tags[i][k] + tmp2;
        }
        break_tags.push_back(tmp);
        break_tags.push_back(tmp2);
      }
    }
  }

  std::sort(break_tags.begin(), break_tags.end());
  auto last = std::unique(break_tags.begin(), break_tags.end());
  break_tags.erase(last, break_tags.end());
  return break_tags;
}

int TagFilterProcessor::compTagScore(const std::string & seq, const std::vector<std::string> & tags) {
  std::vector<int> hit_list(seq.length());
  std::fill(hit_list.begin(), hit_list.end(), 0);
  for (size_t i = 0; i < tags.size(); i++) {
    std::regex rx(tags[i]);
    for (auto it = std::sregex_iterator(seq.begin(), seq.end(), rx);
         it != std::sregex_iterator(); ++it) {
      hit_list[it->position()]++;
    }
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

void TagFilterProcessor::processDB() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();

  LOG_DEBUG("db_file_name " << db_file_name);

  ProteoformPtrVec raw_forms
      = proteoform_factory::readFastaToProteoformPtrVec(db_file_name,
                                                        prsm_para_ptr->getFixModPtrVec());

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

      generateGapEdges(peaks);
      std::vector<std::vector<PeakNodePtr> > componentsFromGraph = getComponentsFromGraph(peaks);
      std::vector<std::string> seq_tag_set = getTags(componentsFromGraph);
      LOG_DEBUG("seq_tag_set.size " << seq_tag_set.size());
      for (size_t i = 0; i < seq_tag_set.size(); i++) {
        LOG_DEBUG("Tag " << i << ": " << seq_tag_set[i]);
      }

      for (size_t k = 0; k < raw_forms.size(); k++) {
        int score = compTagScore(raw_forms[k]->getFastaSeqPtr()->getRawSeq(), seq_tag_set);
        if (score > 0) {
          SimplePrsmPtr simple_prsm = std::make_shared<SimplePrsm>(deconv_ms_ptr_vec[0]->getMsHeaderPtr(),
                                                                   group_spec_num,
                                                                   raw_forms[k], score);
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
