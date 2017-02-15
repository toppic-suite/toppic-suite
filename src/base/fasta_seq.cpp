// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include "base/logger.hpp"
#include "base/ptm_base.hpp"
#include "base/fasta_seq.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

FastaSeq::FastaSeq(const std::string &name_line, 
                   const std::string &ori_seq) {
  int space_pos = name_line.find(" ");
  name_ = name_line.substr(0, space_pos);
  desc_ = name_line.substr(space_pos + 1);
  //rmChar is moved to getAcidPtmPairVec  
  //seq_ = rmChar(ori_seq);
  seq_ = ori_seq;
  compAcidPtmPairVec();
}

FastaSeq::FastaSeq(const std::string &name, 
                   const std::string &desc, 
                   const std::string &ori_seq): 
    name_(name),
    desc_(desc) {
      //seq_ = rmChar(ori_seq);
      seq_ = ori_seq;
      compAcidPtmPairVec();
    }


/** process fasta string and remove unknown letters */
std::string FastaSeq::rmChar(const std::string &ori_seq) {
  std::string seq = "";
  for (size_t i = 0; i < ori_seq.length(); i++) {
    char c = ori_seq.at(i);
    if ((c < 'A' || c > 'Z') ) {
      continue;
    }
    char r = c;
    if (c == 'B') {
      r = 'D';
    } else if (c == 'Z') {
      r = 'E';
    } else if (c == 'X') {
      r = 'A';
    } else if (c == 'J') {
      r = 'I';
    }
    seq = seq + r;
  }
  if (ori_seq != seq) { 
    LOG_INFO( "Reading sequence. Unknown letter occurred. ");
  }
  return seq;
}

void FastaSeq::compAcidPtmPairVec() {
  //LOG_DEBUG("start get acid ptm pair " );
  size_t pos = 0;
  int count = 0;
  while (pos < seq_.length()) {
    std::string acid_one_letter = seq_.substr(pos, 1);
    std::string ptm_str;
    if (pos + 1 >= seq_.length() || seq_.at(pos+1) != '[') {
      ptm_str = PtmBase::getEmptyPtmPtr()->getAbbrName(); 
      pos = pos + 1;
    }
    else {
      //LOG_DEBUG("next letter " << seq_.at(pos+1) << " " << (seq_.at(pos+1) == '['));
      int bracket_pos = seq_.find_first_of("]", pos+1);
      ptm_str = seq_.substr(pos+2, bracket_pos - pos - 2);
      pos = bracket_pos + 1;
    }
    bool valid = true;
    char c = acid_one_letter.at(0);
    if ((c < 'A' || c > 'Z') ) {
      valid = false;
    }
    if (c == 'B') {
      acid_one_letter = "D";
    } else if (c == 'Z') {
      acid_one_letter = "E";
    } else if (c == 'X') {
      acid_one_letter = "A";
    } else if (c == 'J') {
      acid_one_letter = "I";
    }
    if (valid) {
      std::pair<std::string, std::string> pair (acid_one_letter, ptm_str);
      //LOG_DEBUG("count " << count << " acid " << acid_one_letter << " ptm " << ptm_str);
      count++;
      acid_ptm_pair_vec_.push_back(pair);
    }
  }
  //LOG_DEBUG("end get acid ptm pair " );
}

std::string FastaSeq::getString(const std::pair<std::string,std::string> &str_pair) {
  std::string result = str_pair.first;
  std::string ptm_str = str_pair.second;
  if (ptm_str != PtmBase::getEmptyPtmPtr()->getAbbrName()) {
    result = result + "[" + ptm_str + "]";
  }
  return result;
}

std::string FastaSeq::getString(const StringPairVec &str_pair_vec) {
  std::string result;
  for (size_t i = 0; i < str_pair_vec.size(); i++) {
    result = result +getString(str_pair_vec[i]);
  }
  return result;
}

void FastaSeq::appendNameDescToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = FastaSeq::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "seq_name", name_.c_str());
  xml_doc->addElement(element, "seq_desc", desc_.c_str());
  parent->appendChild(element);
}

std::string FastaSeq::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = XmlDomUtil::getChildValue(element, "seq_name", 0);
  return name;
}

std::string FastaSeq::getDescFromXml(xercesc::DOMElement * element) {
  std::string desc = XmlDomUtil::getChildValue(element, "seq_desc", 0);
  return desc;
}

}
