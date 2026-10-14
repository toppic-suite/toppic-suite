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

#include "seq/proteoform.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "common/base/mod_base.hpp"
#include "common/base/mod_type.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/ptm_base.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "seq/alter_type.hpp"
#include "seq/proteoform_factory.hpp"
#include "seq/residue_seq.hpp"

namespace toppic {

Proteoform::Proteoform(const FastaSeqPtr& fasta_seq_ptr,
                       const ProtModPtr& prot_mod_ptr, int start_pos,
                       int end_pos, const ResSeqPtr& res_seq_ptr,
                       const MassShiftPtrVec& mass_shift_ptr_vec)
    : fasta_seq_ptr_(fasta_seq_ptr),
      prot_mod_ptr_(prot_mod_ptr),
      trunc_n_term_mod_ptr_(ModBase::getNTermNoneModPtr()),
      start_pos_(start_pos),
      end_pos_(end_pos),
      residue_seq_ptr_(res_seq_ptr),
      proteo_cluster_id_(-1),
      prot_id_(-1),
      prot_cluster_id_(-1),
      mass_shift_list_(mass_shift_ptr_vec) {
  bp_spec_ptr_ = std::make_shared<BpSpec>(res_seq_ptr);
  std::sort(mass_shift_list_.begin(), mass_shift_list_.end(),
            MassShift::cmpPosInc);
}

Proteoform::Proteoform(XmlDOMElement element,
                       const FastaIndexReaderPtr& reader_ptr,
                       const ModPtrVec& fix_mod_list) {
  std::string seq_element_name = FastaSeq::getXmlElementName();
  XmlDOMElement seq_element =
      xml_dom_util::getChildElement(element, seq_element_name.c_str(), 0);
  std::string seq_name = FastaSeq::getNameFromXml(seq_element);
  std::string seq_desc = FastaSeq::getDescFromXml(seq_element);

  ProteoformPtr form_ptr = proteoform_factory::readFastaToProteoformPtr(
      reader_ptr, seq_name, seq_desc, fix_mod_list);
  parseXml(element, form_ptr);
}

void Proteoform::parseXml(XmlDOMElement element,
                          const ProteoformPtr& form_ptr) {
  start_pos_ = xml_dom_util::getIntChildValue(element, "start_pos", 0);
  end_pos_ = xml_dom_util::getIntChildValue(element, "end_pos", 0);
  proteo_cluster_id_ =
      xml_dom_util::getIntChildValue(element, "proteo_cluster_id", 0);
  proteo_inte_ = xml_dom_util::getDoubleChildValue(element, "proteo_inte", 0);
  prot_id_ = xml_dom_util::getIntChildValue(element, "prot_id", 0);
  prot_cluster_id_ =
      xml_dom_util::getIntChildValue(element, "prot_cluster_id", 0);

  // Get protein N-terminal modification
  std::string pm_element_name = ProtMod::getXmlElementName();
  XmlDOMElement pm_element =
      xml_dom_util::getChildElement(element, pm_element_name.c_str(), 0);
  prot_mod_ptr_ = ProtModBase::getProtModPtrFromXml(pm_element);

  // Add N-terminal modification
  fasta_seq_ptr_ = form_ptr->getFastaSeqPtr();
  residue_seq_ptr_ =
      form_ptr->getResSeqPtr()->getSubResidueSeq(start_pos_, end_pos_);

  ModPtr n_term_mod_ptr = prot_mod_ptr_->getModPtr();
  if (!ModBase::isNTermNoneModPtr(n_term_mod_ptr)) {
    // check only the amino acid to determine whether the modification can be
    // added
    if (residue_seq_ptr_->getLen() >= 1 &&
        ModBase::isNTermNoneModPtr(residue_seq_ptr_->getNModPtr()) &&
        n_term_mod_ptr->getOriResiduePtr()->getAminoAcidPtr() ==
            residue_seq_ptr_->getResiduePtr(0)->getAminoAcidPtr()) {
      residue_seq_ptr_->setNModPtr(n_term_mod_ptr);
    }
  }

  bp_spec_ptr_ = std::make_shared<BpSpec>(residue_seq_ptr_);

  // Get the truncation N-terminal modification (files written before it was
  // added have none)
  trunc_n_term_mod_ptr_ = ModBase::getNTermNoneModPtr();
  if (xml_dom_util::getChildCount(element, "trunc_n_term_mod") > 0) {
    XmlDOMElement trunc_element =
        xml_dom_util::getChildElement(element, "trunc_n_term_mod", 0);
    std::string mod_element_name = Mod::getXmlElementName();
    XmlDOMElement mod_element = xml_dom_util::getChildElement(
        trunc_element, mod_element_name.c_str(), 0);
    if (!setTruncNTermModPtr(ModBase::getModPtrFromXml(mod_element))) {
      LOG_ERROR("The truncation N-terminal modification of a proteoform of "
                << fasta_seq_ptr_->getName() << " cannot be added!");
    }
  }

  // the classification, in the PrSM files of the report stage only
  if (xml_dom_util::getChildCount(element, "classification") > 0) {
    classification_ =
        xml_dom_util::getChildValue(element, "classification", 0);
  }

  // Parse mass shifts
  std::string shift_name = MassShift::getXmlElementName();
  std::string shift_list_name = shift_name + "_list";
  XmlDOMElement list_element =
      xml_dom_util::getChildElement(element, shift_list_name.c_str(), 0);
  int len = xml_dom_util::getChildCount(list_element, shift_name.c_str());

  for (int i = 0; i < len; i++) {
    XmlDOMElement shift_element =
        xml_dom_util::getChildElement(list_element, shift_name.c_str(), i);
    mass_shift_list_.push_back(std::make_shared<MassShift>(shift_element));
  }

  // the alternative mass shift lists (files written before they were
  // added have none)
  std::string alt_list_name = getAlternativeXmlElementName();
  if (xml_dom_util::getChildCount(element, alt_list_name.c_str()) > 0) {
    XmlDOMElement alt_element =
        xml_dom_util::getChildElement(element, alt_list_name.c_str(), 0);
    int alt_num =
        xml_dom_util::getChildCount(alt_element, shift_list_name.c_str());
    for (int k = 0; k < alt_num; k++) {
      XmlDOMElement alt_list =
          xml_dom_util::getChildElement(alt_element, shift_list_name.c_str(), k);
      int n = xml_dom_util::getChildCount(alt_list, shift_name.c_str());
      MassShiftPtrVec shifts;
      for (int i = 0; i < n; i++) {
        XmlDOMElement shift_element =
            xml_dom_util::getChildElement(alt_list, shift_name.c_str(), i);
        shifts.push_back(std::make_shared<MassShift>(shift_element));
      }
      alternative_mass_shift_list_.push_back(shifts);
    }
  }
}

// Get the mass of the modified proteoform
double Proteoform::getMass() {
  double mass = getResSeqPtr()->getSeqMass();
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    // only unexpected and variable changes need to to added
    if (mass_shift_list_[i]->getTypePtr() == AlterType::UNEXPECTED ||
        mass_shift_list_[i]->getTypePtr() == AlterType::VARIABLE) {
      mass += mass_shift_list_[i]->getMassShift();
    }
  }
  return mass;
}

bool Proteoform::hasTruncNTermMod() const {
  return !ModBase::isNTermNoneModPtr(trunc_n_term_mod_ptr_);
}

bool Proteoform::setTruncNTermModPtr(const ModPtr& mod_ptr) {
  if (ModBase::isNTermNoneModPtr(mod_ptr)) {
    if (hasTruncNTermMod()) {
      trunc_n_term_mod_ptr_ = mod_ptr;
      residue_seq_ptr_ = std::make_shared<ResidueSeq>(
          residue_seq_ptr_->getResidues(), mod_ptr);
      bp_spec_ptr_ = std::make_shared<BpSpec>(residue_seq_ptr_);
    }
    return true;
  }
  // the first two residues of a protein are those of its N-terminal forms
  if (prot_mod_ptr_ != ProtModBase::getProtModPtr_NONE() || start_pos_ < 2 ||
      residue_seq_ptr_->getLen() == 0 ||
      mod_ptr->getModTypePtr() != ModType::N_TERM ||
      mod_ptr->getOriResiduePtr()->getAminoAcidPtr() !=
          residue_seq_ptr_->getResiduePtr(0)->getAminoAcidPtr()) {
    return false;
  }
  trunc_n_term_mod_ptr_ = mod_ptr;
  // a new residue sequence, as the old one may be shared
  residue_seq_ptr_ =
      std::make_shared<ResidueSeq>(residue_seq_ptr_->getResidues(), mod_ptr);
  bp_spec_ptr_ = std::make_shared<BpSpec>(residue_seq_ptr_);
  return true;
}

double Proteoform::getMinusWaterMass() {
  double mass = getMass();
  return mass - mass_constant::getWaterMass();
}

PtmPtrVec Proteoform::getPtmVec(const AlterTypePtr& type) {
  PtmPtrVec ptm_vec;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    if (mass_shift_list_[i]->getTypePtr() != type) {
      continue;
    }
    const AlterPtrVec& change_vec = mass_shift_list_[i]->getAlterPtrVec();
    for (size_t k = 0; k < change_vec.size(); k++) {
      ModPtr m = change_vec[k]->getModPtr();
      if (m != nullptr) {
        PtmPtr p = m->getModResiduePtr()->getPtmPtr();
        if (p != nullptr && !PtmBase::isEmptyPtmPtr(p)) {
          ptm_vec.push_back(p);
        }
      }
    }
  }
  return ptm_vec;
}

ProteoformTypePtr Proteoform::getProteoformType() {
  int trunc_len = prot_mod_ptr_->getTruncPtr()->getTruncLen();
  bool is_prefix = false;
  if (start_pos_ == trunc_len) {
    is_prefix = true;
  }

  bool is_suffix = false;
  if (end_pos_ == fasta_seq_ptr_->getAcidPtmPairLen() - 1) {
    is_suffix = true;
  }

  if (is_prefix) {
    if (is_suffix) {
      return ProteoformType::COMPLETE;
    } else {
      return ProteoformType::PREFIX;
    }
  } else {
    if (is_suffix) {
      return ProteoformType::SUFFIX;
    } else {
      return ProteoformType::INTERNAL;
    }
  }
}

int Proteoform::getAlterNum(const AlterTypePtr& type_ptr) {
  int n = 0;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    if (mass_shift_list_[i]->getTypePtr() == type_ptr) {
      n = n + mass_shift_list_[i]->getAlterNum();
    }
  }
  return n;
}

int Proteoform::getVarPtmNum() {
  return getAlterNum(AlterType::VARIABLE) + (hasTruncNTermMod() ? 1 : 0);
}

MassShiftPtrVec Proteoform::getMassShiftPtrVec(const AlterTypePtr& type_ptr) {
  MassShiftPtrVec shift_ptr_vec;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    if (mass_shift_list_[i]->getTypePtr() == type_ptr) {
      shift_ptr_vec.push_back(mass_shift_list_[i]);
    }
  }
  return shift_ptr_vec;
}

void Proteoform::addMassShiftPtrVec(const MassShiftPtrVec& new_shift_ptr_vec) {
  mass_shift_list_.insert(mass_shift_list_.end(), new_shift_ptr_vec.begin(),
                          new_shift_ptr_vec.end());
}

// get several segments without unexpected PTMs from a proteoform
SeqSegmentPtrVec Proteoform::getSeqSegmentPtrVec() {
  MassShiftPtrVec shifts;
  double mass_shift_sum = 0;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    AlterTypePtr type_ptr = mass_shift_list_[i]->getTypePtr();
    if (type_ptr == AlterType::UNEXPECTED || type_ptr == AlterType::VARIABLE) {
      shifts.push_back(mass_shift_list_[i]);
      mass_shift_sum += mass_shift_list_[i]->getMassShift();
    }
  }

  SeqSegmentPtrVec segments;
  double n_shift = 0;
  double c_shift = mass_shift_sum;
  int left = 0;
  for (size_t i = 0; i < shifts.size(); i++) {
    int right = shifts[i]->getLeftBpPos();
    SeqSegmentPtr segment_ptr =
        std::make_shared<SeqSegment>(left, right, n_shift, c_shift);
    segments.push_back(segment_ptr);
    left = shifts[i]->getRightBpPos();
    n_shift = n_shift + shifts[i]->getMassShift();
    c_shift = c_shift - shifts[i]->getMassShift();
  }
  int right = residue_seq_ptr_->getLen();
  SeqSegmentPtr segment_ptr =
      std::make_shared<SeqSegment>(left, right, n_shift, c_shift);
  segments.push_back(segment_ptr);
  return segments;
}

namespace {

// Local function used by getProteoformMatchSeq
void updateMatchSeq(const MassShiftPtrVec& shifts,
                    std::vector<std::string>& left_strings,
                    std::vector<std::string>& right_strings) {
  for (size_t i = 0; i < shifts.size(); i++) {
    MassShiftPtr shift = shifts[i];
    int left_pos = shift->getLeftBpPos();
    int right_pos = shift->getRightBpPos();

    // if it is N-terminimal acetylation
    if (left_pos == 0 && right_pos == 1 &&
        shift->getTypePtr() == AlterType::PROTEIN_VARIABLE &&
        shift->getAlterPtr(0)->getModPtr()->getModResiduePtr()->getPtmPtr() ==
            PtmBase::getPtmPtr_Acetylation()) {
      left_strings[left_pos] = "[Acetyl]-" + left_strings[left_pos];
    } else {
      left_strings[left_pos] = "(" + left_strings[left_pos];
      right_strings[right_pos] += ")";
      right_strings[right_pos] =
          right_strings[right_pos] + "[" + shift->getAnnoStr() + "]";
    }
  }
}

}  // namespace

std::string Proteoform::getAlterStr(const AlterTypePtr& type_ptr) {
  MassShiftPtrVec shifts = getMassShiftPtrVec(type_ptr);
  std::string result;
  for (size_t i = 0; i < shifts.size(); i++) {
    MassShiftPtr shift = shifts[i];
    int left_pos = shift->getLeftBpPos();
    int right_pos = shift->getRightBpPos();
    if ((left_pos + 1) == right_pos) {
      result = result + shift->getAnnoStr() + ":[" +
               std::to_string(left_pos + 1) + "]";
    } else {
      result = result + shift->getAnnoStr() + ":[" +
               std::to_string(left_pos + 1) + "-" + std::to_string(right_pos) +
               "]";
    }
    if (i < shifts.size() - 1) {
      result = result + ";";
    }
  }
  return result;
}

std::string Proteoform::getVarPtmStr() {
  std::string result = getAlterStr(AlterType::VARIABLE);
  if (hasTruncNTermMod()) {
    std::string n_term_str =
        trunc_n_term_mod_ptr_->getModResiduePtr()->getPtmPtr()->getAbbrName() +
        ":[N-term]";
    result = result == "" ? n_term_str : n_term_str + ";" + result;
  }
  return result;
}

std::string Proteoform::getProteoformMatchSeq() {
  std::string mid_string = residue_seq_ptr_->toAcidString();
  std::sort(mass_shift_list_.begin(), mass_shift_list_.end(),
            MassShift::cmpPosInc);

  std::vector<std::string> left_strings(mid_string.size() + 1, "");
  std::vector<std::string> right_strings(mid_string.size() + 1, "");

  MassShiftPtrVec input_shifts = getMassShiftPtrVec(AlterType::INPUT);
  updateMatchSeq(input_shifts, left_strings, right_strings);

  MassShiftPtrVec fixed_shifts = getMassShiftPtrVec(AlterType::FIXED);
  updateMatchSeq(fixed_shifts, left_strings, right_strings);

  MassShiftPtrVec protein_var_shifts =
      getMassShiftPtrVec(AlterType::PROTEIN_VARIABLE);
  updateMatchSeq(protein_var_shifts, left_strings, right_strings);

  MassShiftPtrVec var_shifts = getMassShiftPtrVec(AlterType::VARIABLE);
  updateMatchSeq(var_shifts, left_strings, right_strings);

  MassShiftPtrVec unexpected_shifts = getMassShiftPtrVec(AlterType::UNEXPECTED);
  updateMatchSeq(unexpected_shifts, left_strings, right_strings);

  std::string result = "";
  for (size_t i = 0; i < mid_string.length(); i++) {
    result =
        result + right_strings[i] + left_strings[i] + mid_string.substr(i, 1);
  }
  // last break;
  result = result + right_strings[mid_string.length()];

  // the truncation N-terminal modification, before everything on the first
  // residue: [Acetyl]-(K)[Acetyl]... has both an N-terminal and a side chain
  // acetylation on K
  if (hasTruncNTermMod()) {
    PtmPtr ptm_ptr = trunc_n_term_mod_ptr_->getModResiduePtr()->getPtmPtr();
    result = "[" + ptm_ptr->getAbbrName() + "]-" + result;
  }

  return result;
}

std::string Proteoform::getPrevAminoAcid() {
  const StringPairVec& string_pairs = fasta_seq_ptr_->getAcidPtmPairVec();
  std::string prefix = "-";
  if (start_pos_ > 0) {
    prefix = string_pairs[start_pos_ - 1].first;
  }
  return prefix;
}

std::string Proteoform::getNextAminoAcid() {
  const StringPairVec& string_pairs = fasta_seq_ptr_->getAcidPtmPairVec();
  std::string suffix = "-";
  if (end_pos_ < static_cast<int>(string_pairs.size()) - 1) {
    suffix = string_pairs[end_pos_ + 1].first;
  }
  return suffix;
}

void Proteoform::appendXml(XmlDOMDocument* xml_doc, XmlDOMElement parent) {
  std::string element_name = getXmlElementName();
  XmlDOMElement element = xml_doc->addElement(parent, element_name.c_str());
  fasta_seq_ptr_->appendNameDescToXml(xml_doc, element);
  prot_mod_ptr_->appendNameToXml(xml_doc, element);
  std::string str = std::to_string(start_pos_);
  xml_doc->addElement(element, "start_pos", str.c_str());
  str = std::to_string(end_pos_);
  xml_doc->addElement(element, "end_pos", str.c_str());
  str = std::to_string(proteo_cluster_id_);
  xml_doc->addElement(element, "proteo_cluster_id", str.c_str());
  str = str_util::toString(proteo_inte_);
  xml_doc->addElement(element, "proteo_inte", str.c_str());
  str = std::to_string(prot_id_);
  xml_doc->addElement(element, "prot_id", str.c_str());
  str = std::to_string(prot_cluster_id_);
  xml_doc->addElement(element, "prot_cluster_id", str.c_str());
  // The following three elements are used in prsm_str
  str = std::to_string(getVarPtmNum());
  xml_doc->addElement(element, "variable_ptm_num", str.c_str());
  str = std::to_string(getAlterNum(AlterType::UNEXPECTED));
  xml_doc->addElement(element, "unexpected_ptm_num", str.c_str());
  str = getProteoformMatchSeq();
  xml_doc->addElement(element, "proteo_match_seq", str.c_str());
  str = residue_seq_ptr_->toAcidString();
  xml_doc->addElement(element, "proteo_db_seq", str.c_str());
  // written only when there is one, so that the files of proteoforms
  // without one stay the same
  if (hasTruncNTermMod()) {
    XmlDOMElement trunc_element =
        xml_doc->addElement(element, "trunc_n_term_mod");
    trunc_n_term_mod_ptr_->appendToXml(xml_doc, trunc_element);
  }
  // written only when set (the report stage), so that the PrSM files before
  // it stay the same
  if (classification_ != "") {
    xml_doc->addElement(element, "classification", classification_.c_str());
  }

  element_name = MassShift::getXmlElementName() + "_list";
  XmlDOMElement cl = xml_doc->addElement(element, element_name.c_str());
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    mass_shift_list_[i]->appendXml(xml_doc, cl);
  }
  if (!alternative_mass_shift_list_.empty()) {
    std::string alt_list_name = getAlternativeXmlElementName();
    XmlDOMElement alt_element =
        xml_doc->addElement(element, alt_list_name.c_str());
    for (size_t k = 0; k < alternative_mass_shift_list_.size(); k++) {
      XmlDOMElement alt_list =
          xml_doc->addElement(alt_element, element_name.c_str());
      for (size_t i = 0; i < alternative_mass_shift_list_[k].size(); i++) {
        alternative_mass_shift_list_[k][i]->appendXml(xml_doc, alt_list);
      }
    }
  }
}

std::string Proteoform::getAlternativePtmStr() {
  std::string result;
  for (size_t k = 0; k < alternative_mass_shift_list_.size(); k++) {
    MassShiftPtrVec shifts;
    for (size_t i = 0; i < alternative_mass_shift_list_[k].size(); i++) {
      AlterTypePtr type_ptr = alternative_mass_shift_list_[k][i]->getTypePtr();
      if (type_ptr == AlterType::VARIABLE || type_ptr == AlterType::UNEXPECTED) {
        shifts.push_back(alternative_mass_shift_list_[k][i]);
      }
    }
    std::sort(shifts.begin(), shifts.end(), MassShift::cmpPosInc);
    std::string alt;
    for (size_t i = 0; i < shifts.size(); i++) {
      int left_pos = shifts[i]->getLeftBpPos();
      int right_pos = shifts[i]->getRightBpPos();
      if (i > 0) alt += ";";
      alt += shifts[i]->getAnnoStr() + ":[" + std::to_string(left_pos + 1);
      if (left_pos + 1 != right_pos) alt += "-" + std::to_string(right_pos);
      alt += "]";
    }
    if (k > 0) result += " | ";
    result += alt;
  }
  return result;
}

std::string Proteoform::getMIScore() {
  std::string mi_score = "";

  const StringPairVec& string_pairs = fasta_seq_ptr_->getAcidPtmPairVec();

  MassShiftPtrVec mass_shift_vec = getMassShiftPtrVec(AlterType::VARIABLE);
  for (size_t i = 0; i < mass_shift_vec.size(); i++) {
    if (mass_shift_vec[i]->getAlterPtr(0)->getLocalAnno() == nullptr) continue;

    const std::vector<double>& scr_vec =
        mass_shift_vec[i]->getAlterPtr(0)->getLocalAnno()->getScrVec();
    int left_db_bp = mass_shift_vec[i]->getLeftBpPos() + start_pos_;
    int right_db_bp = mass_shift_vec[i]->getRightBpPos() + start_pos_;
    mi_score = mi_score +
               mass_shift_vec[i]
                   ->getAlterPtr(0)
                   ->getLocalAnno()
                   ->getPtmPtr()
                   ->getAbbrName() +
               "[";

    for (int j = left_db_bp; j < right_db_bp; j++) {
      std::string acid_letter = string_pairs[j].first;
      double scr = std::floor(scr_vec[j - left_db_bp] * 1000) / 10;
      if (scr == 100) scr = 99.9;
      if (scr == 0) continue;

      mi_score =
          mi_score + acid_letter + std::to_string(j - start_pos_ + 1) + ":";
      std::stringstream ss;
      ss << std::fixed << std::setprecision(1) << scr;
      mi_score = mi_score + ss.str() + "%";
      if (j != right_db_bp - 1) {
        mi_score += ";";
      }
    }
    mi_score += "]";

    if (i != mass_shift_vec.size() - 1) {
      mi_score += ";";
    }
  }

  if (mi_score == "") {
    mi_score = "-";
  }
  return mi_score;
}

std::string Proteoform::getClassification(int prot_hit_num,
                                          const ModPtrVec& var_mod_ptrs) {
  bool multiple_prots = prot_hit_num > 1;
  if (getAlterNum(AlterType::UNEXPECTED) > 0) {
    return multiple_prots ? "5" : "4";
  }
  if (hasUndecidedVarPtmSite(var_mod_ptrs)) {
    return multiple_prots ? "3" : "2A";
  }
  return multiple_prots ? "2D" : "1";
}

bool Proteoform::hasUndecidedVarPtmSite(const ModPtrVec& var_mod_ptrs) {
  MassShiftPtrVec shift_ptrs = getMassShiftPtrVec(AlterType::VARIABLE);
  for (size_t i = 0; i < shift_ptrs.size(); i++) {
    int left = shift_ptrs[i]->getLeftBpPos();
    int right = shift_ptrs[i]->getRightBpPos();
    if (right - left <= 1) {
      continue;
    }
    // the PTMs of the shift and their numbers of alterations
    PtmPtrVec ptm_ptrs;
    std::vector<int> alter_nums;
    const AlterPtrVec& alter_ptrs = shift_ptrs[i]->getAlterPtrVec();
    for (size_t j = 0; j < alter_ptrs.size(); j++) {
      if (alter_ptrs[j]->getModPtr() == nullptr) {
        return true;
      }
      const PtmPtr& ptm_ptr =
          alter_ptrs[j]->getModPtr()->getModResiduePtr()->getPtmPtr();
      size_t k = 0;
      while (k < ptm_ptrs.size() && !ptm_ptrs[k]->isSame(ptm_ptr)) {
        k++;
      }
      if (k == ptm_ptrs.size()) {
        ptm_ptrs.push_back(ptm_ptr);
        alter_nums.push_back(0);
      }
      alter_nums[k]++;
    }
    // the candidate residues of each PTM, none a candidate of two
    std::vector<bool> used(right - left, false);
    for (size_t k = 0; k < ptm_ptrs.size(); k++) {
      int cand_num = 0;
      for (int pos = left; pos < right; pos++) {
        const ResiduePtr& res_ptr = residue_seq_ptr_->getResiduePtr(pos);
        bool modifiable = false;
        for (size_t m = 0; m < var_mod_ptrs.size(); m++) {
          if (var_mod_ptrs[m]->getModResiduePtr()->getPtmPtr()->isSame(
                  ptm_ptrs[k]) &&
              res_ptr->isSame(var_mod_ptrs[m]->getOriResiduePtr())) {
            modifiable = true;
            break;
          }
        }
        if (!modifiable) {
          continue;
        }
        if (used[pos - left]) {
          return true;
        }
        used[pos - left] = true;
        cand_num++;
      }
      if (cand_num != alter_nums[k]) {
        return true;
      }
    }
  }
  return false;
}

}  // namespace toppic
