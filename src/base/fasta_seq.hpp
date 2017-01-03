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


#ifndef PROT_BASE_FASTA_SEQ_HPP_
#define PROT_BASE_FASTA_SEQ_HPP_

#include <memory>
#include <string>
#include <vector>

#include "base/residue.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {


class FastaSeq {
 public:
  FastaSeq(const std::string &name_line, const std::string &ori_seq);

  FastaSeq(const std::string &name, const std::string &desc, 
           const std::string &ori_seq);

  std::string getName() {return name_;}

  std::string getDesc() {return desc_;}

  std::string getRawSeq() {return seq_;}

  StringPairVec getAcidPtmPairVec() {return acid_ptm_pair_vec_;}

  int getAcidPtmPairLen() {return acid_ptm_pair_vec_.size();}

  //int getLen() {return seq_.length();}

  static std::string getXmlElementName() {return "fasta_seq";}

  void appendNameDescToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getNameFromXml(xercesc::DOMElement* element);

  static std::string getDescFromXml(xercesc::DOMElement* element);

  static std::string getString(const std::pair<std::string,std::string> &str_pair);

  static std::string getString(const StringPairVec &str_pair_vec);

  static std::string rmChar(const std::string &ori_seq);

 private:
  std::string name_;
  std::string desc_;
  std::string seq_;
  StringPairVec acid_ptm_pair_vec_;

  void compAcidPtmPairVec();
}; 

typedef std::shared_ptr<FastaSeq> FastaSeqPtr;

}

#endif
