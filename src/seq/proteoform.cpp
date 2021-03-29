//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/mod_base.hpp"
#include "common/base/prot_mod_base.hpp"
#include "seq/alter_type.hpp"
#include "seq/fasta_index_reader.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

Proteoform::Proteoform(FastaSeqPtr fasta_seq_ptr,
                       ProtModPtr prot_mod_ptr, 
                       int start_pos, int end_pos,
                       ResSeqPtr res_seq_ptr, 
                       const MassShiftPtrVec & mass_shift_ptr_vec):
    fasta_seq_ptr_(fasta_seq_ptr),
    prot_mod_ptr_(prot_mod_ptr),
    start_pos_(start_pos),
    end_pos_(end_pos),
    residue_seq_ptr_(res_seq_ptr),
    proteo_cluster_id_(-1),
    prot_id_(-1),
    mass_shift_list_(mass_shift_ptr_vec) {
      bp_spec_ptr_ = std::make_shared<BpSpec>(res_seq_ptr);
      std::sort(mass_shift_list_.begin(), mass_shift_list_.end(), 
                MassShift::cmpPosInc);
    }

Proteoform::Proteoform(XmlDOMElement* element, FastaIndexReaderPtr reader_ptr,
                       const ModPtrVec &fix_mod_list) {
  std::string seq_element_name = FastaSeq::getXmlElementName();
  XmlDOMElement* seq_element 
      = xml_dom_util::getChildElement(element, seq_element_name.c_str(), 0);
  std::string seq_name = FastaSeq::getNameFromXml(seq_element);
  std::string seq_desc = FastaSeq::getDescFromXml(seq_element);

  ProteoformPtr form_ptr 
      = proteoform_factory::readFastaToProteoformPtr(reader_ptr, seq_name,
                                                     seq_desc, fix_mod_list);
  parseXml(element, form_ptr);
}

void Proteoform::parseXml(XmlDOMElement* element, ProteoformPtr form_ptr) {
  start_pos_ = xml_dom_util::getIntChildValue(element, "start_pos", 0);
  end_pos_ = xml_dom_util::getIntChildValue(element, "end_pos", 0);
  proteo_cluster_id_ = xml_dom_util::getIntChildValue(element, "proteo_cluster_id", 0);
  prot_id_ = xml_dom_util::getIntChildValue(element, "prot_id", 0);
  variable_ptm_num_ = xml_dom_util::getIntChildValue(element, "variable_ptm_num", 0);

  // Get protein N-terminal modification
  std::string pm_element_name = ProtMod::getXmlElementName();
  XmlDOMElement* pm_element 
      = xml_dom_util::getChildElement(element, pm_element_name.c_str(), 0);
  prot_mod_ptr_ = ProtModBase::getProtModPtrFromXml(pm_element);

  // Add N-terminal modification
  fasta_seq_ptr_ = form_ptr->getFastaSeqPtr();
  residue_seq_ptr_ = form_ptr->getResSeqPtr()->getSubResidueSeq(start_pos_, end_pos_);

  ModPtr mod_ptr = prot_mod_ptr_->getModPtr();
  if (!ModBase::isNoneModPtr(mod_ptr)) {
    if (residue_seq_ptr_->getLen() >= 1 
        && mod_ptr->getOriResiduePtr() == residue_seq_ptr_->getResiduePtr(0)) {
      ResiduePtrVec residues = residue_seq_ptr_->getResidues();
      residues[0] = mod_ptr->getModResiduePtr();
      residue_seq_ptr_ = std::make_shared<ResidueSeq>(residues);
    }
  }

  bp_spec_ptr_ = std::make_shared<BpSpec>(residue_seq_ptr_);

  // Parse mass shifts
  std::string shift_name = MassShift::getXmlElementName();
  std::string shift_list_name = shift_name + "_list";
  XmlDOMElement* list_element 
      = xml_dom_util::getChildElement(element, shift_list_name.c_str(), 0);
  int len = xml_dom_util::getChildCount(list_element, shift_name.c_str());

  for (int i = 0; i < len; i++) {
    XmlDOMElement* shift_element 
        = xml_dom_util::getChildElement(list_element, shift_name.c_str(), i);
    mass_shift_list_.push_back(std::make_shared<MassShift>(shift_element));
  }
}

// Get the mass of the modified proteoform
double Proteoform::getMass() {
  double mass = getResSeqPtr()->getSeqMass();
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    // only unexpected and variable changes need to to added
    if (mass_shift_list_[i]->getTypePtr() == AlterType::UNEXPECTED
        || mass_shift_list_[i]->getTypePtr() == AlterType::VARIABLE) {
      mass += mass_shift_list_[i]->getMassShift();
    }
  }
  return mass;
}

PtmPtrVec Proteoform::getPtmVec(AlterTypePtr type) {
  PtmPtrVec ptm_vec;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    if (mass_shift_list_[i]->getTypePtr() != type) {
      continue;
    }
    AlterPtrVec change_vec = mass_shift_list_[i]->getAlterPtrVec();
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

int Proteoform::getMassShiftNum(AlterTypePtr type_ptr) {
  int n = 0;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    if (mass_shift_list_[i]->getTypePtr() == type_ptr) {
      n++;
    }
  }
  return n;
}

MassShiftPtrVec Proteoform::getMassShiftPtrVec(AlterTypePtr type_ptr) {
  MassShiftPtrVec shift_ptr_vec;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    if (mass_shift_list_[i]->getTypePtr() == type_ptr) {
      shift_ptr_vec.push_back(mass_shift_list_[i]);
    }
  }
  return shift_ptr_vec;
}

void Proteoform::addMassShiftPtrVec(const MassShiftPtrVec & new_shift_ptr_vec) {
  mass_shift_list_.insert(mass_shift_list_.end(), 
                          new_shift_ptr_vec.begin(), 
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
    SeqSegmentPtr segment_ptr = std::make_shared<SeqSegment>(left, right, n_shift, c_shift);
    segments.push_back(segment_ptr);
    left = shifts[i]->getRightBpPos();
    n_shift = n_shift + shifts[i]->getMassShift();
    c_shift = c_shift - shifts[i]->getMassShift();
  }
  int right = residue_seq_ptr_->getLen();
  SeqSegmentPtr segment_ptr = std::make_shared<SeqSegment>(left, right, n_shift, c_shift);
  segments.push_back(segment_ptr);
  return segments;
}

// Local function used by getProteinMatchSeq
void updateMatchSeq(const MassShiftPtrVec & shifts,
                    std::vector<std::string> &left_strings,
                    std::vector<std::string> &right_strings) {
  for (size_t i = 0; i < shifts.size(); i++) {
    int left_pos = shifts[i]->getLeftBpPos();
    left_strings[left_pos] = "(" + left_strings[left_pos];

    int right_pos = shifts[i]->getRightBpPos();
    right_strings[right_pos] +=  ")";
    right_strings[right_pos] = right_strings[right_pos] 
        + "[" + shifts[i]->getAnnoStr() + "]";
  }
}

std::string Proteoform::getProteinMatchSeq() {
  StringPairVec string_pairs = fasta_seq_ptr_->getAcidPtmPairVec();
  std::string mid_string = residue_seq_ptr_->toAcidString();
  std::sort(mass_shift_list_.begin(), mass_shift_list_.end(), MassShift::cmpPosInc);

  std::vector<std::string> left_strings(mid_string.size() + 1, "");
  std::vector<std::string> right_strings(mid_string.size() + 1, "");

  MassShiftPtrVec input_shifts = getMassShiftPtrVec(AlterType::INPUT);
  updateMatchSeq(input_shifts, left_strings, right_strings);

  MassShiftPtrVec fixed_shifts = getMassShiftPtrVec(AlterType::FIXED);
  updateMatchSeq(fixed_shifts, left_strings, right_strings);

  MassShiftPtrVec protein_var_shifts = getMassShiftPtrVec(AlterType::PROTEIN_VARIABLE);
  updateMatchSeq(protein_var_shifts, left_strings, right_strings);

  MassShiftPtrVec var_shifts = getMassShiftPtrVec(AlterType::VARIABLE);
  updateMatchSeq(var_shifts, left_strings, right_strings);

  MassShiftPtrVec unexpected_shifts = getMassShiftPtrVec(AlterType::UNEXPECTED);
  updateMatchSeq(unexpected_shifts, left_strings, right_strings);

  std::string result = "";
  for (size_t i = 0; i < mid_string.length(); i++) {
    result = result + right_strings[i] + left_strings[i] + mid_string.substr(i, 1);
  }
  // last break;
  result = result + right_strings[mid_string.length()];

  std::string prefix = "";
  if (start_pos_ > 0) {
    prefix = string_pairs[start_pos_-1].first;
  }
  std::string suffix = "";
  if (end_pos_ < static_cast<int>(string_pairs.size()) - 1) {
    suffix = string_pairs[end_pos_+1].first;
  }

  return prefix + "." + result + "." + suffix;
}

void Proteoform::appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  fasta_seq_ptr_->appendNameDescToXml(xml_doc, element);
  prot_mod_ptr_->appendNameToXml(xml_doc, element);
  std::string str = str_util::toString(start_pos_);
  xml_doc->addElement(element, "start_pos", str.c_str());
  str = str_util::toString(end_pos_);
  xml_doc->addElement(element, "end_pos", str.c_str());
  str = str_util::toString(proteo_cluster_id_);
  xml_doc->addElement(element, "proteo_cluster_id", str.c_str());
  str = str_util::toString(prot_id_);
  xml_doc->addElement(element, "prot_id", str.c_str());
  str = str_util::toString(variable_ptm_num_);
  xml_doc->addElement(element, "variable_ptm_num", str.c_str());
  str = str_util::toString(getMassShiftNum(AlterType::UNEXPECTED));
  xml_doc->addElement(element, "unexpected_ptm_num", str.c_str());

  element_name = MassShift::getXmlElementName() + "_list";
  XmlDOMElement* cl = xml_doc->createElement(element_name.c_str());
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    mass_shift_list_[i]->appendXml(xml_doc, cl);
  }
  element->appendChild(cl);
  parent->appendChild(element);
}

std::string Proteoform::getMIScore() {
  std::string mi_score = "";

  StringPairVec string_pairs = fasta_seq_ptr_->getAcidPtmPairVec();

  MassShiftPtrVec mass_shift_vec = getMassShiftPtrVec(AlterType::UNEXPECTED);
  for (size_t i = 0; i < mass_shift_vec.size(); i++) {
    if (mass_shift_vec[i]->getAlterPtr(0)->getLocalAnno() == nullptr)
      continue;

    std::vector<double> scr_vec 
        = mass_shift_vec[i]->getAlterPtr(0)->getLocalAnno()->getScrVec();
    int left_db_bp = mass_shift_vec[i]->getLeftBpPos() + start_pos_;
    int right_db_bp = mass_shift_vec[i]->getRightBpPos() + start_pos_;
    mi_score = mi_score 
        + mass_shift_vec[i]->getAlterPtr(0)->getLocalAnno()->getPtmPtr()->getAbbrName() 
        + "[";

    for (int j = left_db_bp; j < right_db_bp; j++) {
      std::string acid_letter = string_pairs[j].first;
      double scr = std::floor(scr_vec[j - left_db_bp] * 1000) / 10;
      if (scr == 100) scr = 99.9;
      if (scr == 0) continue;

      mi_score = mi_score + acid_letter + str_util::toString(j + 1) + ":";
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

}  // namespace toppic

