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


#include <vector>
#include <algorithm>
#include <string>

#include "base/logger.hpp"
#include "base/proteoform_factory.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "spec/ms.hpp"
#include "spec/extend_ms_factory.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm.hpp"

namespace prot {

Prsm::Prsm(ProteoformPtr proteoform_ptr, const DeconvMsPtrVec &deconv_ms_ptr_vec,
           double adjusted_prec_mass, SpParaPtr sp_para_ptr):
    adjusted_prec_mass_(adjusted_prec_mass),
    proteoform_ptr_(proteoform_ptr),
    deconv_ms_ptr_vec_(deconv_ms_ptr_vec) {
      MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getMsHeaderPtr();
      spectrum_id_ = header_ptr->getId();
      spectrum_scan_ = header_ptr->getScansString();
      precursor_id_ = header_ptr->getPrecId();
      prec_feature_id_ = header_ptr->getFeatureId();
      prec_feature_inte_ = header_ptr->getFeatureInte();
      spectrum_num_ = deconv_ms_ptr_vec.size();
      ori_prec_mass_ = header_ptr->getPrecMonoMass();
      init(sp_para_ptr);
    }

Prsm::Prsm(xercesc::DOMElement* element, FastaIndexReaderPtr reader_ptr,
           const ModPtrVec &fix_mod_list) {
  parseXml(element);
  std::string form_elem_name = Proteoform::getXmlElementName();
  xercesc::DOMElement* form_element
      = xml_dom_util::getChildElement(element, form_elem_name.c_str(), 0);
  proteoform_ptr_ = std::make_shared<Proteoform>(form_element, reader_ptr, fix_mod_list);
}

Prsm::Prsm(const Prsm &obj) {
  adjusted_prec_mass_ = obj.adjusted_prec_mass_;
  proteoform_ptr_ = obj.proteoform_ptr_;
  deconv_ms_ptr_vec_ = obj.deconv_ms_ptr_vec_;
  file_name_ = obj.file_name_;
  spectrum_id_ = obj.spectrum_id_;
  spectrum_scan_ = obj.spectrum_scan_;
  precursor_id_ = obj.precursor_id_;
  prec_feature_id_ = obj.prec_feature_id_;
  prec_feature_inte_ = obj.prec_feature_inte_;
  spectrum_num_ = obj.spectrum_num_;
  ori_prec_mass_ = obj.ori_prec_mass_;
  match_peak_num_ = obj.match_peak_num_;
  match_fragment_num_ = obj.match_fragment_num_;
}

void Prsm::init(SpParaPtr sp_para_ptr) {
  refine_ms_three_vec_
      = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec_, sp_para_ptr, adjusted_prec_mass_);
  initScores(sp_para_ptr);
}

void Prsm::initMatchNum(double min_mass) {
  PeakIonPairPtrVec pairs =
      peak_ion_pair_util::genePeakIonPairs(proteoform_ptr_, refine_ms_three_vec_, min_mass);

  match_peak_num_ = 0;
  match_fragment_num_ = 0;
  TheoPeakPtr prev_ion(nullptr);
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getTheoPeakPtr() != prev_ion) {
      prev_ion = pairs[i]->getTheoPeakPtr();
      match_fragment_num_ += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
  std::sort(pairs.begin(), pairs.end(), PeakIonPair::cmpRealPeakPosInc);
  DeconvPeakPtr prev_deconv_peak(nullptr);
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getRealPeakPtr()->getBasePeakPtr() != prev_deconv_peak) {
      prev_deconv_peak = pairs[i]->getRealPeakPtr()->getBasePeakPtr();
      match_peak_num_ += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
}

void Prsm::initScores(SpParaPtr sp_para_ptr) {
  match_fragment_num_ = 0;
  match_peak_num_ = 0;
  for (size_t i = 0; i < refine_ms_three_vec_.size(); i++) {
    // refined one
    PeakIonPairPtrVec pairs =
        peak_ion_pair_util::genePeakIonPairs(proteoform_ptr_, refine_ms_three_vec_[i],
                                             sp_para_ptr->getMinMass());
    match_fragment_num_ += peak_ion_pair_util::compMatchFragNum(pairs);
    match_peak_num_ += peak_ion_pair_util::compMatchPeakNum(pairs);
  }
}

xercesc::DOMElement* Prsm::toXmlElement(XmlDOMDocument* xml_doc) {
  std::string element_name = Prsm::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "file_name", file_name_.c_str());
  std::string str = string_util::convertToString(prsm_id_);
  xml_doc->addElement(element, "prsm_id", str.c_str());
  str = string_util::convertToString(spectrum_id_);
  xml_doc->addElement(element, "spectrum_id", str.c_str());
  xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
  str = string_util::convertToString(precursor_id_);
  xml_doc->addElement(element, "precursor_id", str.c_str());
  str = string_util::convertToString(prec_feature_id_);
  xml_doc->addElement(element, "precursor_feature_id", str.c_str());
  str = string_util::convertToString(prec_feature_inte_);
  xml_doc->addElement(element, "precursor_feature_inte", str.c_str());
  str = string_util::convertToString(spectrum_num_);
  xml_doc->addElement(element, "spectrum_number", str.c_str());
  str = string_util::convertToString(ori_prec_mass_);
  xml_doc->addElement(element, "ori_prec_mass", str.c_str());
  str = string_util::convertToString(adjusted_prec_mass_);
  xml_doc->addElement(element, "adjusted_prec_mass", str.c_str());
  str = string_util::convertToString(fdr_);
  xml_doc->addElement(element, "fdr", str.c_str());
  str = string_util::convertToString(proteoform_fdr_);
  xml_doc->addElement(element, "proteoform_fdr", str.c_str());
  str = string_util::convertToString(match_peak_num_);
  xml_doc->addElement(element, "match_peak_num", str.c_str());
  str = string_util::convertToString(match_fragment_num_);
  xml_doc->addElement(element, "match_fragment_num", str.c_str());
  str = string_util::convertToString(getNormMatchFragNum());
  xml_doc->addElement(element, "norm_match_fragment_num", str.c_str());
  proteoform_ptr_->appendXml(xml_doc, element);
  if (extreme_value_ptr_ != nullptr) {
    extreme_value_ptr_->appendXml(xml_doc, element);
  }
  return element;
}

void Prsm::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = toXmlElement(xml_doc);
  parent->appendChild(element);
}

void Prsm::parseXml(xercesc::DOMElement *element) {
  file_name_ = xml_dom_util::getChildValue(element, "file_name", 0);
  prsm_id_ = xml_dom_util::getIntChildValue(element, "prsm_id", 0);
  spectrum_id_ = xml_dom_util::getIntChildValue(element, "spectrum_id", 0);
  spectrum_scan_ = xml_dom_util::getChildValue(element, "spectrum_scan", 0);
  precursor_id_ = xml_dom_util::getIntChildValue(element, "precursor_id", 0);
  prec_feature_id_ = xml_dom_util::getIntChildValue(element, "precursor_feature_id", 0);
  prec_feature_inte_ = xml_dom_util::getDoubleChildValue(element, "precursor_feature_inte", 0);
  spectrum_num_ = xml_dom_util::getIntChildValue(element, "spectrum_number", 0);
  ori_prec_mass_ = xml_dom_util::getDoubleChildValue(element, "ori_prec_mass", 0);
  adjusted_prec_mass_ = xml_dom_util::getDoubleChildValue(element, "adjusted_prec_mass", 0);
  fdr_ = xml_dom_util::getDoubleChildValue(element, "fdr", 0);
  proteoform_fdr_ = xml_dom_util::getDoubleChildValue(element, "proteoform_fdr", 0);
  match_peak_num_ = xml_dom_util::getDoubleChildValue(element, "match_peak_num", 0);
  match_fragment_num_ = xml_dom_util::getDoubleChildValue(element, "match_fragment_num", 0);

  int prob_count = xml_dom_util::getChildCount(element, "extreme_value");
  if (prob_count != 0) {
    xercesc::DOMElement* prob_element
        = xml_dom_util::getChildElement(element, "extreme_value", 0);
    extreme_value_ptr_ = std::make_shared<ExtremeValue>(prob_element);
  }
}

double Prsm::getEValue() {
  if (extreme_value_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  } else {
    return extreme_value_ptr_->getEValue();
  }
}

double Prsm::getPValue() {
  if (extreme_value_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  } else {
    return extreme_value_ptr_->getPValue();
  }
}

double Prsm::getOneProtProb() {
  if (extreme_value_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  } else {
    return extreme_value_ptr_->getOneProtProb();
  }
}

/* this function is tempory for testing mass graph alignment */
double Prsm::getNormMatchFragNum() {
  int var_change_num = proteoform_ptr_->getVariablePtmNum();
  int unexp_change_num = proteoform_ptr_->getMassShiftNum(MassShiftType::UNEXPECTED);
  int start_pos = proteoform_ptr_->getStartPos();
  int end_pos = proteoform_ptr_->getEndPos();
  double score = match_fragment_num_ - 2 * var_change_num - 4 * unexp_change_num;
  if (start_pos == 0 || start_pos == 1) {
    score = score + 2;
  }
  if (end_pos == static_cast<int>(proteoform_ptr_->getFastaSeqPtr()->getRawSeq().length()) - 1) {
    score = score + 2;
  }
  return score;
}

// sort by the number of matched fragments, then the number of matched peaks
bool Prsm::cmpMatchFragmentDecMatchPeakDec(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getMatchFragNum() > b->getMatchFragNum()) {
    return true;
  } else if (a->getMatchFragNum() < b->getMatchFragNum()) {
    return false;
  }
  return a->getMatchPeakNum() > b->getMatchPeakNum();
}

// sort by number of matched fragment ions, then start position
bool Prsm::cmpMatchFragDecStartPosInc(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getMatchFragNum() > b->getMatchFragNum()) {
    return true;
  } else if (a->getMatchFragNum() == b->getMatchFragNum()) {
    return a->getProteoformPtr()->getStartPos() < b->getProteoformPtr()->getStartPos();
  }
  return false;
}

// sort by the order of spectrum id, the precursor id
bool Prsm::cmpSpectrumIdIncPrecursorIdInc(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getSpectrumId() < b->getSpectrumId()) {
    return true;
  } else if (a->getSpectrumId() > b->getSpectrumId()) {
    return false;
  } else {
    if (a->getPrecursorId() < b->getPrecursorId()) {
      return true;
    }
    return false;
  }
}

// sort by spectrum id then match ions
bool Prsm::cmpSpectrumIdIncMatchFragDec(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getSpectrumId() < b->getSpectrumId()) {
    return true;
  } else if (a->getSpectrumId() > b->getSpectrumId()) {
    return false;
  } else {
    if (a->getMatchFragNum() > b->getMatchFragNum()) {
      return true;
    }
    return false;
  }
}

// sort by spectrum id then evalue
bool Prsm::cmpSpectrumIdIncEvalueInc(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getSpectrumId() < b->getSpectrumId()) {
    return true;
  } else if (a->getSpectrumId() > b->getSpectrumId()) {
    return false;
  } else {
    if (a->getEValue() < b->getEValue()) {
      return true;
    }
    return false;
  }
}

}  // namespace prot
