
#include <locale> 
#include "tagfinder.hpp"
#include "peak_node.hpp"

namespace prot{

Tagfinder::Tagfinder(PrsmParaPtr prsm_para_ptr, int gap){
  gap_ = gap;
  prsm_para_ptr_ = prsm_para_ptr;
  ResiduePtrVec residue_list = 
      ResidueUtil::convertStrToResiduePtrVec("ARNDCEQGHILKMFPSTWYVUO", 
                                             prsm_para_ptr->getFixModPtrVec());
  for (size_t i = 0; i < residue_list.size(); i++) {
    mass_list_.push_back(std::make_pair(residue_list[i]->getAcidPtr()->getOneLetter(), residue_list[i]->getMass()));
  }

  if (gap_ > 1) {
    for (size_t i = 0; i < residue_list.size(); i++) {
      for (size_t j = i + 1; j < residue_list.size(); j++) {
        mass_list_.push_back(std::make_pair(
                std::to_string(residue_list[i]->getMass() + residue_list[j]->getMass()), 
                residue_list[i]->getMass() + residue_list[j]->getMass()));
        if (gap_ > 2) {
          for (size_t k = j + 1; k < residue_list.size(); k++) {
            mass_list_.push_back(std::make_pair(
                    std::to_string(residue_list[i]->getMass() + residue_list[j]->getMass() + residue_list[k]->getMass()), 
                    residue_list[i]->getMass() + residue_list[j]->getMass() + residue_list[k]->getMass()));
          }
        }

      }
    }
  }
}

void Tagfinder::process(){

  std::string spec_file_name = prsm_para_ptr_->getSpectrumFileName();

  MsAlignReader sp_reader(spec_file_name, 1);
  SpectrumSetPtr spec_set_ptr;
  SpParaPtr sp_para_ptr = prsm_para_ptr_->getSpParaPtr();
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    if(spec_set_ptr->isValid()){
      DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
      PrmMsPtrVec ms_six_ptr_vec = spec_set_ptr->getMsSixPtrVec();
      PeakTolerancePtr tole_ptr = prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
      PrmPeakPtrVec prm_peak_vec = PrmMs::getPrmPeakPtrs(ms_six_ptr_vec, tole_ptr);
      std::vector<PeakNodePtr> peaks = generateGapEdges(prm_peak_vec, gap_);
      std::vector<std::vector<PeakNodePtr> > componentsFromGraph = getComponentsFromGraph(peaks);
      std::vector<std::string> tags = getTags(componentsFromGraph);
      for (size_t i = 0; i < tags.size(); i++) {
        std::cout << spec_set_ptr->getSpecId() << " " << tags[i] << std::endl;
      }
    }
  }
}

std::vector<PeakNodePtr> Tagfinder::generateGapEdges(const PrmPeakPtrVec & prm_peaks, int gap) {
  std::vector<PeakNodePtr> peaks(prm_peaks.size());
  for (size_t i = 0; i < peaks.size(); i++) {
    peaks[i] = std::make_shared<PeakNode>(prm_peaks[i]);
    peaks[i]->clearEdges();
  }

  std::sort(peaks.begin(), peaks.end(), 
            [](PeakNodePtr a, PeakNodePtr b){return a->getMonoMass() < b->getMonoMass();});

  for (size_t i = 0; i < peaks.size() - 1; i++) {
    for (size_t j = i + 1; j < peaks.size(); j++) {
      double err = peaks[i]->getMonoMass() * prsm_para_ptr_->getErrorTolerance() / 1000000.0;
      double diff = peaks[j]->getMonoMass() - peaks[i]->getMonoMass();
      for (size_t k = 0; k < mass_list_.size(); k++) {
        if (std::abs(diff - mass_list_[k].second) < err) {
          peaks[i]->addNext(peaks[j]);
          break;
        }
      }
    }
  }
  return peaks;
}

std::vector<std::vector<PeakNodePtr>> getComponentsFromGraph(std::vector<PeakNodePtr> peaks) {
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
    int componentId = p->getComponentId();
    if (!componentDone[componentId]) {
      std::vector<PeakNodePtr> component;
      for (size_t j = 0; j < peaks.size(); j++) {
        if (peaks[j]->getComponentId() == componentId) {
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

std::vector<std::string> Tagfinder::getTags(std::vector<std::vector<PeakNodePtr>> componentsFromGraph) {
  std::vector<std::string> tags;
  std::locale loc;
  for (size_t i = 0; i < componentsFromGraph.size(); i++) {
    std::vector<PeakNodePtr> tagPeaks = findBestTag(componentsFromGraph[i]);
    std::string superTag = "";
    if (tagPeaks.size() < TAG_SIZE + 1) {
      continue;
    }
    for (size_t j = 0; j < tagPeaks.size() - 1; j++) {
      double err = tagPeaks[j]->getMonoMass() * prsm_para_ptr_->getErrorTolerance() / 1000000.0;
      double diff = tagPeaks[j + 1]->getMonoMass() - tagPeaks[j]->getMonoMass();
      for (size_t k = 0; k < mass_list_.size(); k++) {
        if (std::abs(diff - mass_list_[k].second) < err) {
          if (std::isalpha(mass_list_[k].first[0],loc)) {
            superTag += mass_list_[k].first;
          } else {
            superTag += "[" + mass_list_[k].first + "]";
          }
          break;
        }
      }
    }
    tags.push_back(superTag);
  }
  return tags;
}

std::vector<PeakNodePtr> findBestTag(PeakNodePtr peak, std::vector<PeakNodePtr> best, int len, std::vector<PeakNodePtr> & prefix){
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

std::vector<PeakNodePtr> findBestTag(std::vector<PeakNodePtr> peaks){
  std::vector<PeakNodePtr> bestTag;

  for (size_t i = 0; i < peaks.size(); i++) {
    peaks[i]->setMaxPrefix(-1);
  }

  std::vector<PeakNodePtr> prefix;
  for (int i = 0; i < 500; i++) {
    prefix.push_back(nullptr);
  }

  for (size_t i = 0; i < peaks.size(); i++) {
    bestTag = findBestTag(peaks[i], bestTag, 0, prefix);
  }
  return bestTag;
}

}
