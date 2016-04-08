#include "base/acid_base.hpp"
#include "base/residue_util.hpp"
#include "base/mod_util.hpp"
#include "spec/theo_peak_factory.hpp"
#include "tag_filter.hpp"

namespace prot{

double getPeptideMass(const std::string & seq) {
  double m = 0;
  for (size_t i = 0; i < seq.length(); i++) {
    m += AcidBase::getAcidPtrByOneLetter(seq.substr(i, 1))->getMonoMass();
  }
  return m;
}

SeqTag TagFilter::geneSeqTag(ProteoformPtr proteoform) {
  SeqTag seq_tag;
  std::string seq = proteoform->getFastaSeqPtr()->getSeq(); 
  for (size_t i = 0; i < seq.length() - 1; i++) {
    std::string s = seq.substr(i, 2);
    std::vector<double> mass;
    for (size_t j = 0; j <= i; j++) {
      mass.push_back(getPeptideMass(seq.substr(0, j)));
    }
    seq_tag[s].push_back(std::make_pair(i, mass));
  }

  return seq_tag;
}

TagFilter::TagFilter(const ProteoformPtrVec &proteo_ptrs,
                     TagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    seq_tag_vec_.push_back(geneSeqTag(proteo_ptrs[i]));
    seq_name_vec_.push_back(proteo_ptrs[i]->getFastaSeqPtr()->getName());
  }

  ResiduePtrVec residue_list = 
      ResidueUtil::convertStrToResiduePtrVec("ARNDCEQGHILKMFPSTWYVUO", 
                                             mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  ModPtrVec var_mod_ptr_vec = ModUtil::readModTxt(mng_ptr_->var_mod_file_name_)[2];

  for (size_t j = 0; j < var_mod_ptr_vec.size(); j++) {
    residue_list.push_back(var_mod_ptr_vec[j]->getModResiduePtr());
  }

  for (size_t i = 0; i < residue_list.size(); i++) {
    residue_mass_list_.push_back(std::make_pair(residue_list[i]->getAcidPtr()->getOneLetter(), residue_list[i]->getMass()));
    for (size_t j = i + 1; j < residue_list.size(); j++) {
      residue_mass_list_.push_back(
          std::make_pair(residue_list[i]->getAcidPtr()->getOneLetter()
                         + residue_list[j]->getAcidPtr()->getOneLetter(),
                         residue_list[i]->getMass() + residue_list[j]->getMass())); 
    }
  }
}

void count(const std::vector<std::pair<int, std::vector<double> > > & seq_tag,
           double min, double max, std::map<int, int> & counter) {
  if (seq_tag.size() == 0) return; 
  for (size_t i = 0; i < seq_tag.size(); i++) {
    for (size_t j = 0; j < seq_tag[i].second.size(); j++) {
      if (seq_tag[i].second[j] >= min && seq_tag[i].second[j] <= max)
        counter[j]++;
    }
  }
}

int compute(SeqTag seq_tag, const std::vector<SpecTag> & spec_tag) {
  std::map<int, int> counter;
  for (size_t i = 0; i < spec_tag.size(); i++) {
    SpecTag t = spec_tag[i];
    std::string tag_seg = t.getSeq();
    count(seq_tag[tag_seg], t.getMinMass(), t.getMaxMass(), counter);
    if (!t.isOrdered()) {
      tag_seg = tag_seg.substr(1, 1) + tag_seg.substr(0, 1);
      count(seq_tag[tag_seg], t.getMinMass(), t.getMaxMass(), counter);
    }
  }

  int max = 0;
  for(std::map<int, int>::iterator it = counter.begin(); it!= counter.end(); it++) {
    if(it->second > max) max = it->second; 
  }
  return max;
}

std::vector<std::string> TagFilter::getBestMatch(const PrmMsPtrVec &ms_ptr_vec) {
  std::vector<std::string> res;
  PeakTolerancePtr tole_ptr = 
      mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  PrmPeakPtrVec prm_peak_vec = PrmMs::getPrmPeakPtrs(ms_ptr_vec, tole_ptr);
  std::vector<SpecTag> spec_tag_vec = geneSpecTag(prm_peak_vec);
  std::vector<std::pair<std::string, int> > counter;
  for (size_t i = 0; i < seq_tag_vec_.size(); i++) {
    counter.push_back(std::make_pair(seq_name_vec_[i], compute(seq_tag_vec_[i], spec_tag_vec))); 
  }

  std::sort(counter.begin(), counter.end(), 
            [](std::pair<std::string, int> a, std::pair<std::string, int> b){
            return a.second < b.second;
            });
  return res;
}

std::vector<SpecTag> TagFilter::geneSpecTag(const PrmPeakPtrVec & prm_peaks) {
  std::vector<SpecTag> tag_vec;

  std::vector<double> mass_list;
  for (size_t i = 0; i < prm_peaks.size(); i++) {
    mass_list.push_back(prm_peaks[i]->getMonoMass());
  }
  std::sort(mass_list.begin(), mass_list.end());
  for (size_t i = 0; i < mass_list.size() - 2; i++) {
    for (size_t j = i + 1; j < mass_list.size() - 1; j++) {
      double diff = mass_list[j] - mass_list[i];
      double err = mass_list[i] * mng_ptr_->prsm_para_ptr_->getErrorTolerance() / 1000000;
      for (size_t r = 0; r < residue_mass_list_.size(); r++) {
        if (std::abs(diff - residue_mass_list_[r].second) < err) {
          if (residue_mass_list_[r].first.length() == 2) {
            SpecTag t(residue_mass_list_[r].first.substr(0,1),
                      residue_mass_list_[r].first.substr(1,1),
                      mass_list[i], false); 
            tag_vec.push_back(t);
          } else if (residue_mass_list_[r].first.length() == 1) {
            for (size_t k = j + 1; k < mass_list.size(); k++) {
              double diff2 = mass_list[k] - mass_list[j]; 
              double err2 = mass_list[j] * mng_ptr_->prsm_para_ptr_->getErrorTolerance() / 1000000;
              for (size_t r2 = 0; r2 < residue_mass_list_.size(); r2++) {
                if (residue_mass_list_[r2].first.length() > 1) continue;
                if (std::abs(diff2 - residue_mass_list_[r2].second) < err2) {
                  SpecTag t(residue_mass_list_[r].first,
                            residue_mass_list_[r2].first,
                            mass_list[i], true);
                  tag_vec.push_back(t);
                }
              }
            }
          } 
        }
      }
    }
  }

  return tag_vec;
}

}
