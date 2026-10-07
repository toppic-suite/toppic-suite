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

#include "seq/proteoform_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "common/base/mod_base.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/prot_mod_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/residue_util.hpp"
#include "common/util/str_util.hpp"
#include "seq/alter_type.hpp"
#include "seq/local_anno.hpp"
#include "seq/mass_shift.hpp"

namespace toppic {
namespace proteoform_util {

ResFreqPtrVec compNTermResidueFreq(const ProteoformPtrVec& prot_mod_forms) {
  std::vector<double> counts;
  ResiduePtrVec residue_list;
  for (size_t i = 0; i < prot_mod_forms.size(); i++) {
    ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();
    if (seq_ptr->getLen() >= 1) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(0);
      int pos = residue_util::findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found
        counts[pos] = counts[pos] + 1;
      } else {
        residue_list.push_back(res_ptr);
        counts.push_back(1);
      }
    }
  }

  double sum = 0;
  for (size_t i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (size_t i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr = std::make_shared<ResidueFreq>(
        residue_list[i]->getAminoAcidPtr(), residue_list[i]->getPtmPtr(),
        counts[i] / sum);
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}

ResFreqPtrVec compResidueFreq(const ResiduePtrVec& residue_list,
                              const ProteoformPtrVec& prot_mod_forms) {
  std::vector<double> counts(residue_list.size(), 0.0);
  for (size_t i = 0; i < prot_mod_forms.size(); i++) {
    ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();
    for (int j = 0; j < seq_ptr->getLen(); j++) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(j);
      int pos = residue_util::findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found
        counts[pos] = counts[pos] + 1;
      }
    }
  }

  double sum = 0;
  for (size_t i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (size_t i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr = std::make_shared<ResidueFreq>(
        residue_list[i]->getAminoAcidPtr(), residue_list[i]->getPtmPtr(),
        counts[i] / sum);
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}

bool isSameSeqAndMass(const ProteoformPtr& a, const ProteoformPtr& b,
                      double ppo) {
  if (a->getSeqName() != b->getSeqName()) {
    return false;
  }

  if (a->getStartPos() != b->getStartPos()) {
    return false;
  }

  if (a->getEndPos() != b->getEndPos()) {
    return false;
  }

  double thresh = a->getMass() * ppo;

  if (std::abs(a->getMass() - b->getMass()) > thresh) {
    return false;
  }

  return true;
}

bool isStrictCompatiablePtmSpecies(const ProteoformPtr& a,
                                   const ProteoformPtr& b, double ppo) {
  if (!isSameSeqAndMass(a, b, ppo)) {
    return false;
  }

  if (a->getMassShiftNum() != b->getMassShiftNum()) {
    return false;
  }

  double shift_tolerance = a->getResSeqPtr()->getSeqMass() * ppo;
  // sort mass shifts
  MassShiftPtrVec a_shift_vec = a->getMassShiftPtrVec();
  MassShiftPtrVec b_shift_vec = b->getMassShiftPtrVec();
  std::sort(a_shift_vec.begin(), a_shift_vec.end(), MassShift::cmpPosInc);
  std::sort(b_shift_vec.begin(), b_shift_vec.end(), MassShift::cmpPosInc);
  for (int i = 0; i < a->getMassShiftNum(); i++) {
    MassShiftPtr ac = a_shift_vec[i];
    MassShiftPtr bc = b_shift_vec[i];
    if (ac->getRightBpPos() <= bc->getLeftBpPos() ||
        bc->getRightBpPos() <= ac->getLeftBpPos()) {
      return false;
    }
    if (std::abs(ac->getMassShift() - bc->getMassShift()) > shift_tolerance) {
      return false;
    }
  }
  return true;
}

ProteoformPtrVec2D divideProteoIntoBlocks(const ProteoformPtrVec& proteo_ptrs,
                                          int db_block_size) {
  size_t start_idx = 0;
  size_t proteo_idx = 0;
  int block_len = 0;
  ProteoformPtrVec2D proteo_blocks;
  while (proteo_idx < proteo_ptrs.size()) {
    int proteo_len = proteo_ptrs[proteo_idx]->getResSeqPtr()->getLen();
    if (block_len + proteo_len < db_block_size) {
      block_len = block_len + proteo_len;
      proteo_idx++;
    } else {
      size_t end_idx = proteo_idx;
      ProteoformPtrVec proteo_in_block;
      for (size_t i = start_idx; i <= end_idx; i++) {
        proteo_in_block.push_back(proteo_ptrs[i]);
      }
      proteo_blocks.push_back(proteo_in_block);
      start_idx = end_idx + 1;
      proteo_idx = end_idx + 1;
      block_len = 0;
    }
  }
  // last block
  if (start_idx < proteo_ptrs.size()) {
    ProteoformPtrVec proteo_in_block;
    for (size_t i = start_idx; i < proteo_ptrs.size(); i++) {
      proteo_in_block.push_back(proteo_ptrs[i]);
    }
    proteo_blocks.push_back(proteo_in_block);
  }
  return proteo_blocks;
}

std::vector<double> getNTermShift(const ProteoformPtr& db_form_ptr,
                                  const ProtModPtrVec& prot_mod_ptrs) {
  std::vector<double> shifts;
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    ResSeqPtr db_res_seq_ptr = db_form_ptr->getResSeqPtr();
    bool valid = prot_mod_util::allowMod(prot_mod_ptrs[i],
                                         db_res_seq_ptr->getResidues());
    if (valid) {
      shifts.push_back(prot_mod_ptrs[i]->getProtShift());
    }
  }
  // sort shifts in the decreasing order
  std::sort(shifts.begin(), shifts.end(), std::greater<>());
  return shifts;
}

std::vector<double> getNTermAcets(const ProteoformPtr& db_form_ptr,
                                  const ProtModPtrVec& prot_mod_ptrs) {
  std::vector<double> shifts;
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    // check if it is acetylation
    if (prot_mod_ptrs[i]->isAcetylation()) {
      ResSeqPtr db_res_seq_ptr = db_form_ptr->getResSeqPtr();
      bool valid = prot_mod_util::allowMod(prot_mod_ptrs[i],
                                           db_res_seq_ptr->getResidues());
      if (valid) {
        shifts.push_back(prot_mod_ptrs[i]->getProtShift());
      }
    }
  }
  return shifts;
}

std::vector<std::vector<double> > getNTermShift2D(
    const ProteoformPtrVec& db_form_ptr_vec,
    const ProtModPtrVec& prot_mod_ptrs) {
  std::vector<std::vector<double> > shifts_2d;
  for (size_t i = 0; i < db_form_ptr_vec.size(); i++) {
    std::vector<double> shifts =
        getNTermShift(db_form_ptr_vec[i], prot_mod_ptrs);
    shifts_2d.push_back(shifts);
  }
  return shifts_2d;
}

std::vector<std::vector<double> > getNTermAcet2D(
    const ProteoformPtrVec& db_form_ptr_vec,
    const ProtModPtrVec& prot_mod_ptrs) {
  std::vector<std::vector<double> > shifts_2d;
  for (size_t i = 0; i < db_form_ptr_vec.size(); i++) {
    std::vector<double> shifts =
        getNTermAcets(db_form_ptr_vec[i], prot_mod_ptrs);
    shifts_2d.push_back(shifts);
  }
  return shifts_2d;
}

namespace {

// A delta mass in ProForma: sign always written, 4 decimals.
std::string proFormaMass(double mass) {
  std::string str = str_util::fixedToString(mass, 4);
  if (str[0] != '-') {
    str = "+" + str;
  }
  return str;
}

// Whether a name can stand in a ProForma tag: no "|" or "#" (the spec
// requires the accession for such CV names) and paired square brackets
// (Cation:Mg[II]).
bool isProFormaName(const std::string& name) {
  if (name.empty() || name.find_first_of("|#") != std::string::npos) {
    return false;
  }
  int depth = 0;
  for (char c : name) {
    if (c == '[') {
      depth++;
    } else if (c == ']' && --depth < 0) {
      return false;
    }
  }
  return depth == 0;
}

// A PTM in ProForma: with a Unimod id, its Unimod name and accession,
// joined with "|" as synonyms, [Oxidation|UNIMOD:35] (the name from
// ptm_base.xml, which holds every Unimod modification, whatever the
// modification file calls the PTM), or the accession alone when the id has
// no name there; without a Unimod id, the mass (the PTM's, or the
// alteration's) with the PTM's name as information, [+12.3457|INFO:MyPTM].
std::string proFormaPtm(const PtmPtr& ptm_ptr, double mass) {
  int unimod_id = ptm_ptr->getUnimodId();
  if (unimod_id > 0) {
    std::string accession = "UNIMOD:" + std::to_string(unimod_id);
    std::string unimod_name = PtmBase::getUnimodName(unimod_id);
    if (isProFormaName(unimod_name)) {
      return unimod_name + "|" + accession;
    }
    return accession;
  }
  std::string str = proFormaMass(mass);
  const std::string& name = ptm_ptr->getAbbrName();
  return isProFormaName(name) ? str + "|INFO:" + name : str;
}

std::string proFormaPtm(const PtmPtr& ptm_ptr) {
  return proFormaPtm(ptm_ptr, ptm_ptr->getMonoMass());
}

// The modification of an alteration in ProForma: the PTM of the local
// annotation or of the mod, the mass of an unexpected shift.
std::string proFormaAlter(const AlterPtr& alter_ptr) {
  if (alter_ptr->getLocalAnno() != nullptr &&
      alter_ptr->getLocalAnno()->getPtmPtr() != nullptr) {
    return proFormaPtm(alter_ptr->getLocalAnno()->getPtmPtr());
  }
  if (alter_ptr->getTypePtr() == AlterType::UNEXPECTED ||
      alter_ptr->getModPtr() == nullptr) {
    return proFormaMass(alter_ptr->getMass());
  }
  PtmPtr ptm_ptr = alter_ptr->getModPtr()->getModResiduePtr()->getPtmPtr();
  if (ptm_ptr == nullptr) {
    return proFormaMass(alter_ptr->getMass());
  }
  // without a Unimod id, the mass of the alteration (the PTM's, or a
  // residue change)
  return proFormaPtm(ptm_ptr, alter_ptr->getMass());
}

// A tag with a label (#s1(0.900), #g1): the label goes after the tag's
// first part, before an INFO part, which would take it in.
std::string labelTag(const std::string& tag, const std::string& label) {
  size_t pos = tag.find('|');
  if (pos == std::string::npos) {
    return tag + label;
  }
  return tag.substr(0, pos) + label + tag.substr(pos);
}

std::string scoreStr(double score) { return str_util::fixedToString(score, 3); }

// The modifications whose position is a range of residues [left, right).
struct ProFormaRange {
  int left;
  int right;
  std::vector<std::string> tags;
};

}  // namespace

std::string getProFormaStr(const ProteoformPtr& form_ptr) {
  std::string acids = form_ptr->getResSeqPtr()->toAcidString();
  int len = acids.size();
  // the modifications after each residue, the N-terminal ones and those of
  // unknown position
  std::vector<std::vector<std::string>> res_tags(len);
  std::vector<std::string> n_term_tags;
  std::vector<std::string> unknown_tags;
  std::vector<ProFormaRange> ranges;
  int group_num = 0;

  if (form_ptr->hasTruncNTermMod()) {
    n_term_tags.push_back(proFormaPtm(
        form_ptr->getTruncNTermModPtr()->getModResiduePtr()->getPtmPtr()));
  }

  MassShiftPtrVec shifts = form_ptr->getMassShiftPtrVec();
  std::stable_sort(shifts.begin(), shifts.end(), MassShift::cmpPosInc);
  for (const MassShiftPtr& shift_ptr : shifts) {
    for (const AlterPtr& alter_ptr : shift_ptr->getAlterPtrVec()) {
      // the alteration's residues, or its mass shift's when the alteration
      // has none within them (TopMG writes an unexpected shift on one
      // residue with an empty alteration range)
      int left = alter_ptr->getLeftBpPos();
      int right = alter_ptr->getRightBpPos();
      if (left >= right || left < shift_ptr->getLeftBpPos() ||
          right > shift_ptr->getRightBpPos()) {
        left = shift_ptr->getLeftBpPos();
        right = shift_ptr->getRightBpPos();
      }
      std::string tag = proFormaAlter(alter_ptr);
      // the protein N-terminal acetylations
      if (alter_ptr->getTypePtr() == AlterType::PROTEIN_VARIABLE &&
          left == 0 && right == 1 && alter_ptr->getModPtr() != nullptr &&
          alter_ptr->getModPtr()->getModResiduePtr()->getPtmPtr() ==
              PtmBase::getPtmPtr_Acetylation()) {
        n_term_tags.push_back(tag);
        continue;
      }
      if (left < 0 || right > len || left >= right) {
        unknown_tags.push_back(tag);
        continue;
      }
      LocalAnnoPtr anno_ptr = alter_ptr->getLocalAnno();
      if (anno_ptr != nullptr && right - left > 1) {
        // MIScore: the site scores of the residues, indexed by the residue
        // position in the proteoform (or from the left end of the range)
        const std::vector<double>& scr_vec = anno_ptr->getScrVec();
        int offset = static_cast<int>(scr_vec.size()) == len ? 0 : left;
        std::vector<int> sites;
        int best = -1;
        for (int p = left; p < right; p++) {
          int idx = p - offset;
          if (idx < 0 || idx >= static_cast<int>(scr_vec.size()) ||
              scr_vec[idx] <= 0) {
            continue;
          }
          sites.push_back(p);
          if (best < 0 || scr_vec[idx] > scr_vec[best - offset]) {
            best = p;
          }
        }
        if (sites.size() == 1) {
          res_tags[best].push_back(tag);
          continue;
        }
        if (sites.size() > 1) {
          std::string label = "#s" + std::to_string(++group_num);
          for (int p : sites) {
            std::string score = "(" + scoreStr(scr_vec[p - offset]) + ")";
            res_tags[p].push_back(p == best ? labelTag(tag, label + score)
                                            : label + score);
          }
          continue;
        }
      }
      if (right - left == 1) {
        res_tags[left].push_back(tag);
        continue;
      }
      // a range: the modifications of one range are written together
      bool added = false;
      for (ProFormaRange& range : ranges) {
        if (range.left == left && range.right == right) {
          range.tags.push_back(tag);
          added = true;
          break;
        }
      }
      if (!added) {
        ranges.push_back(ProFormaRange{left, right, {tag}});
      }
    }
  }

  // ranges may not overlap: the first of overlapping ranges (from the left,
  // the wider one first) stays a range, the others become sets of positions
  std::stable_sort(ranges.begin(), ranges.end(),
                   [](const ProFormaRange& a, const ProFormaRange& b) {
                     if (a.left != b.left) {
                       return a.left < b.left;
                     }
                     return a.right > b.right;
                   });
  std::vector<int> range_begin(len, -1);
  std::vector<int> range_end(len + 1, -1);
  std::vector<bool> covered(len, false);
  for (size_t r = 0; r < ranges.size(); r++) {
    const ProFormaRange& range = ranges[r];
    bool overlap = false;
    for (int p = range.left; p < range.right; p++) {
      if (covered[p]) {
        overlap = true;
        break;
      }
    }
    if (!overlap) {
      for (int p = range.left; p < range.right; p++) {
        covered[p] = true;
      }
      range_begin[range.left] = r;
      range_end[range.right] = r;
      continue;
    }
    for (const std::string& tag : range.tags) {
      std::string label = "#g" + std::to_string(++group_num);
      res_tags[range.left].push_back(labelTag(tag, label));
      for (int p = range.left + 1; p < range.right; p++) {
        res_tags[p].push_back(label);
      }
    }
  }

  std::string result;
  for (const std::string& tag : unknown_tags) {
    result += "[" + tag + "]";
  }
  if (!unknown_tags.empty()) {
    result += "?";
  }
  for (const std::string& tag : n_term_tags) {
    result += "[" + tag + "]";
  }
  if (!n_term_tags.empty()) {
    result += "-";
  }
  for (int i = 0; i < len; i++) {
    if (range_begin[i] >= 0) {
      result += "(";
    }
    result += acids[i];
    for (const std::string& tag : res_tags[i]) {
      result += "[" + tag + "]";
    }
    if (range_end[i + 1] >= 0) {
      result += ")";
      for (const std::string& tag : ranges[range_end[i + 1]].tags) {
        result += "[" + tag + "]";
      }
    }
  }
  return result;
}

}  // namespace proteoform_util
}  // namespace toppic
