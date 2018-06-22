//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include "base/logger.hpp"
#include "base/mass_shift_type.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod_base.hpp"
#include "base/string_util.hpp"
#include "base/fasta_index_reader.hpp"
#include "base/proteoform.hpp"
#include "base/proteoform_factory.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

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
      std::sort(mass_shift_list_.begin(), mass_shift_list_.end(), MassShift::cmpPosInc);
    }

Proteoform::Proteoform(xercesc::DOMElement* element, FastaIndexReaderPtr reader_ptr,
                       const ModPtrVec &fix_mod_list) {
  std::string seq_element_name = FastaSeq::getXmlElementName();
  xercesc::DOMElement* seq_element = xml_dom_util::getChildElement(element, seq_element_name.c_str(), 0);
  std::string seq_name = FastaSeq::getNameFromXml(seq_element);
  std::string seq_desc = FastaSeq::getDescFromXml(seq_element);

  ProteoformPtr form_ptr = proteoform_factory::readFastaToProteoformPtr(reader_ptr, seq_name,
                                                                        seq_desc, fix_mod_list);
  parseXml(element, form_ptr);
}

void Proteoform::parseXml(xercesc::DOMElement* element, ProteoformPtr form_ptr) {
  // LOG_DEBUG("start parse proteoform");
  start_pos_ = xml_dom_util::getIntChildValue(element, "start_pos", 0);
  end_pos_ = xml_dom_util::getIntChildValue(element, "end_pos", 0);
  proteo_cluster_id_ = xml_dom_util::getIntChildValue(element, "proteo_cluster_id", 0);
  prot_id_ = xml_dom_util::getIntChildValue(element, "prot_id", 0);
  variable_ptm_num_ = xml_dom_util::getIntChildValue(element, "variable_ptm_num", 0);

  // LOG_DEBUG("start parse prot mod");
  std::string pm_element_name = ProtMod::getXmlElementName();
  xercesc::DOMElement* pm_element = xml_dom_util::getChildElement(element, pm_element_name.c_str(), 0);
  prot_mod_ptr_ = ProtModBase::getProtModPtrFromXml(pm_element);

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

  // LOG_DEBUG("start parse changes");
  std::string shift_element_name = MassShift::getXmlElementName();

  xercesc::DOMElement* change_list_element = xml_dom_util::getChildElement(element, "mass_shift_list", 0);
  int shift_len = xml_dom_util::getChildCount(change_list_element, shift_element_name.c_str());

  for (int i = 0; i < shift_len; i++) {
    xercesc::DOMElement* shift_element
        = xml_dom_util::getChildElement(change_list_element, shift_element_name.c_str(), i);
    mass_shift_list_.push_back(std::make_shared<MassShift>(shift_element));
  }
  // LOG_DEBUG("end parse proteoform");
}

// get mass of the modified proteoform
double Proteoform::getMass() {
  double mass = getResSeqPtr()->getSeqMass();
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    // only unexpected and variable changes need to to added
    if (mass_shift_list_[i]->getTypePtr() == MassShiftType::UNEXPECTED
        || mass_shift_list_[i]->getTypePtr() == MassShiftType::VARIABLE) {
      mass += mass_shift_list_[i]->getMassShift();
    }
  }
  return mass;
}

PtmPtrVec Proteoform::getPtmVec() {
  PtmPtrVec ptm_vec;

  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    ChangePtrVec change_vec = mass_shift_list_[i]->getChangePtrVec();

    for (size_t k = 0; k < change_vec.size(); k++) {
      PtmPtr p = change_vec[k]->getModPtr()->getModResiduePtr()->getPtmPtr();
      if (p != nullptr && !PtmBase::isEmptyPtmPtr(p)) {
        ptm_vec.push_back(p);
      }
    }
  }

  return ptm_vec;
}

AlignTypePtr Proteoform::getAlignType() {
  int trunc_len = prot_mod_ptr_->getTruncPtr()->getTruncLen();
  // LOG_DEBUG("seq " << getProteinMatchSeq() << " trunc len " << trunc_len << " start pos " << start_pos_);
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
      return AlignType::COMPLETE;
    } else {
      return AlignType::PREFIX;
    }
  } else {
    if (is_suffix) {
      return AlignType::SUFFIX;
    } else {
      return AlignType::INTERNAL;
    }
  }
}

int Proteoform::getMassShiftNum(MassShiftTypePtr ct_ptr) {
  int n = 0;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    if (mass_shift_list_[i]->getTypePtr() == ct_ptr) {
      n++;
    }
  }
  return n;
}

MassShiftPtrVec Proteoform::getMassShiftPtrVec(MassShiftTypePtr ct_ptr) {
  MassShiftPtrVec shift_ptr_vec;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    if (mass_shift_list_[i]->getTypePtr() == ct_ptr) {
      shift_ptr_vec.push_back(mass_shift_list_[i]);
    }
  }
  return shift_ptr_vec;
}

void Proteoform::addMassShiftPtrVec(const MassShiftPtrVec & new_shift_ptr_vec) {
  mass_shift_list_.insert(mass_shift_list_.end(), new_shift_ptr_vec.begin(), new_shift_ptr_vec.end());
}

// get several segments without unexpected PTMs from a proteoform
SegmentPtrVec Proteoform::getSegmentPtrVec() {
  MassShiftPtrVec shifts;
  double mass_shift_sum = 0;
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    MassShiftTypePtr type_ptr = mass_shift_list_[i]->getTypePtr();
    if (type_ptr == MassShiftType::UNEXPECTED || type_ptr == MassShiftType::VARIABLE) {
      shifts.push_back(mass_shift_list_[i]);
      mass_shift_sum += mass_shift_list_[i]->getMassShift();
    }
  }

  SegmentPtrVec segments;
  double n_shift = 0;
  double c_shift = mass_shift_sum;
  int left = 0;
  for (size_t i = 0; i < shifts.size(); i++) {
    int right = shifts[i]->getLeftBpPos();
    SegmentPtr segment_ptr = std::make_shared<Segment>(left, right, n_shift, c_shift);
    segments.push_back(segment_ptr);
    left = shifts[i]->getRightBpPos();
    n_shift = n_shift + shifts[i]->getMassShift();
    c_shift = c_shift - shifts[i]->getMassShift();
  }
  int right = residue_seq_ptr_->getLen();
  SegmentPtr segment_ptr = std::make_shared<Segment>(left, right, n_shift, c_shift);
  segments.push_back(segment_ptr);
  return segments;
}

void updateMatchSeq(const MassShiftPtrVec & shifts,
                    std::vector<std::string> &left_strings,
                    std::vector<std::string> &right_strings) {
  for (size_t i = 0; i < shifts.size(); i++) {
    int left_pos = shifts[i]->getLeftBpPos();
    left_strings[left_pos] = "(" + left_strings[left_pos];

    int right_pos = shifts[i]->getRightBpPos();
    right_strings[right_pos] +=  ")";
    right_strings[right_pos] = right_strings[right_pos] + "[" + shifts[i]->getSeqStr() + "]";
  }
}

std::string Proteoform::getProteinMatchSeq() {
  StringPairVec string_pairs = fasta_seq_ptr_->getAcidPtmPairVec();
  // LOG_DEBUG("string_pairs length " << string_pairs.size() << " string " << FastaSeq::getString(string_pairs));
  std::string mid_string = residue_seq_ptr_->toAcidString();
  // LOG_DEBUG("mid string lenth " << mid_string.length() << " string " << mid_string);
  std::sort(mass_shift_list_.begin(), mass_shift_list_.end(), MassShift::cmpPosInc);

  std::vector<std::string> left_strings(mid_string.size() + 1, "");
  std::vector<std::string> right_strings(mid_string.size() + 1, "");

  // LOG_DEBUG("mass shift update started");
  MassShiftPtrVec input_shifts = getMassShiftPtrVec(MassShiftType::INPUT);
  updateMatchSeq(input_shifts, left_strings, right_strings);

  MassShiftPtrVec fixed_shifts = getMassShiftPtrVec(MassShiftType::FIXED);
  updateMatchSeq(fixed_shifts, left_strings, right_strings);

  MassShiftPtrVec protein_var_shifts = getMassShiftPtrVec(MassShiftType::PROTEIN_VARIABLE);
  updateMatchSeq(protein_var_shifts, left_strings, right_strings);

  MassShiftPtrVec var_shifts = getMassShiftPtrVec(MassShiftType::VARIABLE);
  updateMatchSeq(var_shifts, left_strings, right_strings);

  MassShiftPtrVec unexpected_shifts = getMassShiftPtrVec(MassShiftType::UNEXPECTED);
  updateMatchSeq(unexpected_shifts, left_strings, right_strings);
  // LOG_DEBUG("mass shift update completed");

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

  // LOG_DEBUG("Prefix " << prefix << " result " << result << " suffix length " << suffix.length() << " suffix " << suffix);
  return prefix + "." + result + "." + suffix;
}

std::string Proteoform::toString() {
  std::stringstream s;
  s << "Begin pos: " << start_pos_ << std::endl;
  s << "End pos: " << end_pos_ << std::endl;
  s << "String: " << residue_seq_ptr_->toString();
  return s.str();
}

void Proteoform::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  fasta_seq_ptr_->appendNameDescToXml(xml_doc, element);
  prot_mod_ptr_->appendNameToXml(xml_doc, element);
  std::string str = string_util::convertToString(start_pos_);
  xml_doc->addElement(element, "start_pos", str.c_str());
  str = string_util::convertToString(end_pos_);
  xml_doc->addElement(element, "end_pos", str.c_str());
  str = string_util::convertToString(proteo_cluster_id_);
  xml_doc->addElement(element, "proteo_cluster_id", str.c_str());
  str = string_util::convertToString(prot_id_);
  xml_doc->addElement(element, "prot_id", str.c_str());
  str = string_util::convertToString(variable_ptm_num_);
  xml_doc->addElement(element, "variable_ptm_num", str.c_str());
  str = string_util::convertToString(getMassShiftNum(MassShiftType::UNEXPECTED));
  xml_doc->addElement(element, "unexpected_ptm_num", str.c_str());

  element_name = MassShift::getXmlElementName() + "_list";
  xercesc::DOMElement* cl = xml_doc->createElement(element_name.c_str());
  for (size_t i = 0; i < mass_shift_list_.size(); i++) {
    mass_shift_list_[i]->appendXml(xml_doc, cl);
  }
  element->appendChild(cl);
  parent->appendChild(element);
}

std::string Proteoform::getMIScore() {
  mi_score_ = "";

  StringPairVec string_pairs = fasta_seq_ptr_->getAcidPtmPairVec();

  MassShiftPtrVec mass_shift_vec = getMassShiftPtrVec(MassShiftType::UNEXPECTED);

  for (size_t i = 0; i < mass_shift_vec.size(); i++) {
    if (mass_shift_vec[i]->getChangePtr(0)->getLocalAnno() == nullptr)
      continue;

    std::vector<double> scr_vec = mass_shift_vec[i]->getChangePtr(0)->getLocalAnno()->getScrVec();
    int left_db_bp = mass_shift_vec[i]->getLeftBpPos() + start_pos_;
    int right_db_bp = mass_shift_vec[i]->getRightBpPos() + start_pos_;
    mi_score_ = mi_score_ + mass_shift_vec[i]->getChangePtr(0)->getLocalAnno()->getPtmPtr()->getAbbrName() + "[";

    for (int j = left_db_bp; j < right_db_bp; j++) {
      std::string acid_letter = string_pairs[j].first;
      double scr = std::floor(scr_vec[j - left_db_bp] * 1000) / 10;
      if (scr == 100) scr = 99.9;
      if (scr == 0) continue;

      mi_score_ = mi_score_ + acid_letter + std::to_string(j + 1) + ":";
      std::stringstream ss;
      ss << std::fixed << std::setprecision(1) << scr;
      mi_score_ = mi_score_ + ss.str() + "%";
      if (j != right_db_bp - 1) {
        mi_score_ += ";";
      }
    }
    mi_score_ += "]";

    if (i != mass_shift_vec.size() - 1) {
      mi_score_ += ";";
    }
  }

  if (mi_score_ == "") {
    mi_score_ = "-";
  }
  return mi_score_;
}

}  // namespace prot

