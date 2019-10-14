//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_util.hpp"
#include "seq/fasta_seq.hpp"

namespace toppic {

FastaSeq::FastaSeq(const std::string &name_line,
                   const std::string &ori_seq) {
  int space_pos = name_line.find(" ");
  name_ = name_line.substr(0, space_pos);
  desc_ = name_line.substr(space_pos + 1);
  seq_ = ori_seq;
  compAcidPtmPairVec();
}

FastaSeq::FastaSeq(const std::string &name,
                   const std::string &desc,
                   const std::string &ori_seq):
    name_(name),
    desc_(desc),
    seq_(ori_seq) {
      compAcidPtmPairVec();
    }

// This function need to be tested
FastaSeq::FastaSeq(FastaSeqPtr seq_ptr, int start, int len) {
  name_ = seq_ptr->getName();
  desc_ = seq_ptr->getDesc();
  StringPairVec str_pair_vec = seq_ptr->getAcidPtmPairVec();
  acid_ptm_pair_vec_.insert(acid_ptm_pair_vec_.begin(), 
                            str_pair_vec.begin() + start, 
                            str_pair_vec.begin() + start + len); 
  compRawSeq();
}

void FastaSeq::compAcidPtmPairVec() {
  size_t pos = 0;
  int count = 0;
  while (pos < seq_.length()) {
    std::string acid_one_letter = seq_.substr(pos, 1);
    std::string ptm_str;
    if (pos + 1 >= seq_.length() || seq_.at(pos+1) != '[') {
      ptm_str = PtmBase::getEmptyPtmPtr()->getAbbrName();
      pos = pos + 1;
    } else {
      int bracket_pos = seq_.find_first_of("]", pos+1);
      ptm_str = seq_.substr(pos+2, bracket_pos - pos - 2);
      pos = bracket_pos + 1;
    }
    char c = acid_one_letter.at(0);
    if (residue_util::isValidResidue(c)) {
      acid_one_letter = residue_util::replaceResidueLetter(c);
      std::pair<std::string, std::string> pair(acid_one_letter, ptm_str);
      count++;
      acid_ptm_pair_vec_.push_back(pair);
    }
    else {
      LOG_WARN("The residue " << acid_one_letter << " is invalid!");
    }
  }
}

void FastaSeq::compRawSeq() {
  for (size_t i = 0; i < acid_ptm_pair_vec_.size(); i++) {
    std::string amino_str = acid_ptm_pair_vec_[i].first;
    seq_ = seq_ + amino_str; 
    std::string ptm_str = acid_ptm_pair_vec_[i].second;
    if (ptm_str != PtmBase::getEmptyPtmPtr()->getAbbrName()) {
      seq_ = seq_ + "[" + ptm_str + "]";
    }
  }
}

void FastaSeq::appendNameDescToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = FastaSeq::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "seq_name", name_.c_str());
  xml_doc->addElement(element, "seq_desc", desc_.c_str());
  parent->appendChild(element);
}

std::string FastaSeq::getNameFromXml(XmlDOMElement * element) {
  std::string name = xml_dom_util::getChildValue(element, "seq_name", 0);
  return name;
}

std::string FastaSeq::getDescFromXml(XmlDOMElement * element) {
  std::string desc = xml_dom_util::getChildValue(element, "seq_desc", 0);
  return desc;
}

}  // namespace toppic
