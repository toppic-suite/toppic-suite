//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_COMMON_SEQ_FASTA_SEQ_HPP_
#define TOPPIC_COMMON_SEQ_FASTA_SEQ_HPP_

#include <memory>
#include <string>
#include <vector>

#include "common/seq/amino_acid_replace.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_element.hpp"

namespace toppic {

class XmlDOMDocument;

class FastaSeq;

using FastaSeqPtr = std::shared_ptr<FastaSeq>;

class FastaSeq {
 public:
  FastaSeq(const std::string &name_line, const std::string &ori_seq);

  FastaSeq(const std::string &name, const std::string &desc,
           const std::string &ori_seq);

  FastaSeq(FastaSeqPtr seq_ptr, int start, int len);

  const std::string& getName() const {return name_;}

  const std::string& getDesc() const {return desc_;}

  const std::string& getRawSeq() const {return seq_;}

  std::string getSubSeq(int start, int end) const;

  const StringPairVec& getAcidPtmPairVec() const {return acid_ptm_pair_vec_;}

  std::string getAcidReplaceStr(int bgn, int end) const;

  int getAcidPtmPairLen() const {return acid_ptm_pair_vec_.size();}

  void appendNameDescToXml(XmlDOMDocument* xml_doc, XmlDOMElement parent) const;

  static std::string getXmlElementName() {return "fasta_seq";}

  static std::string getNameFromXml(XmlDOMElement element);

  static std::string getDescFromXml(XmlDOMElement element);

 private:
  std::string name_;

  std::string desc_;

  std::string seq_;

  StringPairVec acid_ptm_pair_vec_;

  AminoAcidReplacePtrVec acid_replace_ptr_vec_;

  void compAcidPtmPairVec();
};

using FastaSeqPtrVec = std::vector<FastaSeqPtr>;

}  // namespace toppic

#endif
