// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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

#include <boost/regex.hpp>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/db_block.hpp"
#include "base/file_util.hpp"
#include "base/proteoform_factory.hpp"
#include "base/residue_util.hpp"
#include "base/mod_util.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_table_writer.hpp"
#include "spec/msalign_util.hpp"
#include "tag_filter_processor.hpp"
#include "peak_node.hpp"
#include "sequence_tag_graph.hpp"

namespace prot{

std::string geneRegex(std::string seq) {
  std::sort(seq.begin(), seq.end());
  std::vector<std::string> perm_list;
  do {
    perm_list.push_back(seq);
  } while(std::next_permutation(seq.begin(), seq.end()));

  std::string rx = "[";
  for (size_t i = 0; i < perm_list.size(); i++) {
    rx = rx + "(" + perm_list[i] + ")";
    if (i != perm_list.size() - 1) {
      rx = rx + "|";
    }
  }
  rx += "]";
  return rx;
}

void TagFilterProcessor::process(){
  ModPtrVec fix_mod = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();
  ResiduePtrVec residue_list = 
      ResidueUtil::convertStrToResiduePtrVec("ARNDCEQGHILKMFPSTWYVUO", fix_mod);
  if (mng_ptr_->residueModFileName_ != "") {
    ModPtrVec var_mod = ModUtil::readModTxt(mng_ptr_->residueModFileName_)[2];
    for (size_t j = 0; j < var_mod.size(); j++) {
      residue_list.push_back(var_mod[j]->getModResiduePtr());
    }
  }

  // generate mass list with variable ptms
  for (size_t i = 0; i < residue_list.size(); i++) {
    mass_list_.push_back(std::make_pair(residue_list[i]->getMass(),
                                        residue_list[i]->getAcidPtr()->getOneLetter()));
    if (mng_ptr_->gap_ > 1) {
      for (size_t j = i; j < residue_list.size(); j++) {
        if (residue_list[i]->getAcidPtr()->getOneLetter() ==
            residue_list[j]->getAcidPtr()->getOneLetter()) {
          mass_list_.push_back(std::make_pair(residue_list[i]->getMass() + residue_list[j]->getMass(),
                                              residue_list[i]->getAcidPtr()->getOneLetter() +
                                              residue_list[j]->getAcidPtr()->getOneLetter()));
        } else {
          mass_list_.push_back(std::make_pair(residue_list[i]->getMass() + residue_list[j]->getMass(),
                                              geneRegex(residue_list[i]->getAcidPtr()->getOneLetter() +
                                                        residue_list[j]->getAcidPtr()->getOneLetter())));
        }
        if (mng_ptr_->gap_ > 2) {
          for (size_t k = j; k < residue_list.size(); k++) {
            if (residue_list[i]->getAcidPtr()->getOneLetter() == 
                residue_list[j]->getAcidPtr()->getOneLetter() &&
                residue_list[k]->getAcidPtr()->getOneLetter() == 
                residue_list[j]->getAcidPtr()->getOneLetter()) {
              mass_list_.push_back(std::make_pair(residue_list[i]->getMass() + residue_list[j]->getMass() +
                                                  residue_list[k]->getMass(),
                                                  residue_list[i]->getAcidPtr()->getOneLetter() +
                                                  residue_list[j]->getAcidPtr()->getOneLetter() +
                                                  residue_list[k]->getAcidPtr()->getOneLetter()));

            } else {
              mass_list_.push_back(std::make_pair(residue_list[i]->getMass() + residue_list[j]->getMass() +
                                                  residue_list[k]->getMass(),
                                                  geneRegex(residue_list[i]->getAcidPtr()->getOneLetter() +
                                                            residue_list[j]->getAcidPtr()->getOneLetter() +
                                                            residue_list[k]->getAcidPtr()->getOneLetter())));
            }
          }
        }
      } 
    }
  }

  std::sort(mass_list_.begin(), mass_list_.end(),
            [](const std::pair<double, std::string> & a,
               const std::pair<double, std::string> & b)->bool
            {
            return a.first < b.first;
            });
  for (size_t i = 0; i < mass_list_.size(); i++) {
    std::cout << mass_list_[i].first << " " << mass_list_[i].second << std::endl;
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
      if (peaks[i]->updateComponentId()){
        done = false;
      }
    }
  } while(!done);

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

  if (len >= (int)best.size()) {
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
    }if (tags[i].size() == mng_ptr_->tag_min_len_) {
      std::string tmp = "";
      std::string tmp2 = "";
      //for (auto const& s : tags[i]) { tmp += s; }
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
        for (size_t k = j; k < j + mng_ptr_->tag_min_len_;k++) {
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
    boost::regex rx(tags[i]);
    for (auto it = boost::sregex_iterator(seq.begin(), seq.end(), rx);
         it != boost::sregex_iterator(); ++it) {
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
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();

  ProteoformPtrVec raw_forms 
      = ProteoformFactory::readFastaToProteoformPtrVec(db_file_name, 
                                                       prsm_para_ptr->getFixModPtrVec());

  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr =  mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(),
                       group_spec_num,
                       sp_para_ptr->getActivationPtr());

  std::string output_file_name = FileUtil::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_;

  SimplePrsmXmlWriter writer(output_file_name);

  SpectrumSetPtr spec_set_ptr;

  while((spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)) != nullptr){
    // basic spectrum graph
    DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
    LOG_DEBUG("Spec ID " <<  deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getId());
    if (spec_set_ptr->isValid()) {
      std::vector<PeakNodePtr> peaks;
      for (size_t k = 0; k < deconv_ms_ptr_vec[0]->size(); k++) {
        peaks.push_back(std::make_shared<PeakNode>(deconv_ms_ptr_vec[0]->getPeakPtr(k)->getMonoMass()));
      }

      std::sort(peaks.begin(), peaks.end(),[](PeakNodePtr a, PeakNodePtr b)->bool
                {
                return a->getMass() < b->getMass();
                });
#ifdef MS_TAG
      // tag generation from Ms-align+Tag finder
      generateGapEdges(peaks);
      std::vector<std::vector<PeakNodePtr> > componentsFromGraph = getComponentsFromGraph(peaks);
      std::vector<std::string> seq_tag_set = getTags(componentsFromGraph);
      std::cout << __FILE__ << " seq_tag_set.size " << seq_tag_set.size() << std::endl;
      for (size_t i = 0; i < seq_tag_set.size(); i++) {
        std::cout << __FILE__ << "Tag" << i << ": " << seq_tag_set[i] << std::endl; 
      }
#endif
#ifdef MS_PATHFINDER
      // tag generation from mspathfinder
      std::shared_ptr<SequenceTagGraph> seq_tag_graph = std::make_shared<SequenceTagGraph>();
      seq_tag_graph->setErrTolerance(mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo());
      seq_tag_graph->SetAminoAcidsArray(mass_list_);
      seq_tag_graph->SetPeakList(peaks);
      seq_tag_graph->SetNodeCount(peaks.size());
      seq_tag_graph->CollectSequenceTagGraphEdges();
      std::vector<std::string> seq_tag_set = seq_tag_graph->FindSequenceTags();
#endif
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
  SimplePrsmReader sim_prsm_reader(output_file_name);
  SimplePrsmPtr sim_prsm_ptr = sim_prsm_reader.readOnePrsm();
  SimplePrsmPtrVec sim_prsm_list;
  while (sim_prsm_ptr != nullptr) {
    sim_prsm_list.push_back(sim_prsm_ptr);
    sim_prsm_ptr = sim_prsm_reader.readOnePrsm();
  }
  sim_prsm_reader.close();
  std::sort(sim_prsm_list.begin(), sim_prsm_list.end(), SimplePrsm::cmpScoreDec);
  SimplePrsmXmlWriter sim_prsm_writer(output_file_name);
  for (size_t i = 0; i < mng_ptr_->top_num_ && i < sim_prsm_list.size(); i++) {
    sim_prsm_writer.write(sim_prsm_list[i]);
  } 
  sim_prsm_writer.close();
#ifdef MS_TAG
  SimplePrsmTableWriterPtr table_out = 
      std::make_shared<SimplePrsmTableWriter>(mng_ptr_->prsm_para_ptr_,
                                              mng_ptr_->output_file_ext_, "TAG_FILTER_TABLE");
#endif
#ifdef MS_PATHFINDER
  SimplePrsmTableWriterPtr table_out = 
      std::make_shared<SimplePrsmTableWriter>(mng_ptr_->prsm_para_ptr_,
                                              mng_ptr_->output_file_ext_, "MS_PATHFINDER_TABLE");
#endif

  table_out->write();
  table_out = nullptr;
}

}
