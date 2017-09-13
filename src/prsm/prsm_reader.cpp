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

#include <string>
#include <vector>

#include "xercesc/framework/MemBufInputSource.hpp"
#include "htslib/faidx.h"

#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "prsm/prsm_reader.hpp"

namespace prot {

PrsmReader::PrsmReader(const std::string &file_name) {
  input_.open(file_name.c_str(), std::ios::in);
}

std::vector<std::string> PrsmReader::readOnePrsmLines() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    line = StringUtil::trim(line);
    // LOG_DEBUG("line " << line);
    if (line ==  "<prsm>") {
      line_list.push_back(line);
    } else if (line == "</prsm>") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    } else if (line == "") {
      continue;
    } else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

PrsmStrPtr PrsmReader::readOnePrsmStr() {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  if (prsm_str_vec.size() == 0) {
    return PrsmStrPtr(nullptr);
  }
  return std::make_shared<PrsmStr>(prsm_str_vec);
}

PrsmPtr PrsmReader::readOnePrsm(FastaIndexReaderPtr reader_ptr,
                                const ModPtrVec fix_mod_list) {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  // LOG_DEBUG("prsm str vec size " << prsm_str_vec.size());
  if (prsm_str_vec.size() == 0) {
    return PrsmPtr(nullptr);
  }
  std::string prsm_str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
  for (size_t i = 0; i < prsm_str_vec.size(); i++) {
    prsm_str += prsm_str_vec[i];
  }
  // LOG_DEBUG("prsm str " << prsm_str);
  xercesc::MemBufInputSource prsm_buf(
      (const XMLByte*)prsm_str.c_str(), prsm_str.size(), "prsm_str (in memory)");

  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  PrsmPtr ptr;
  if (parser) {
    XmlDOMDocument doc(parser, prsm_buf);
    xercesc::DOMElement* root = doc.getDocumentElement();
    ptr = std::make_shared<Prsm>(root, reader_ptr, fix_mod_list);
  }
  // LOG_DEBUG("simple prsm spectrum id " << ptr->getSpectrumId() << " seq name " << ptr->getSeqName());
  return ptr;
}

void PrsmReader::close() {
  input_.close();
}

PrsmStrPtrVec PrsmReader::readAllPrsmStrs(const std::string &input_file_name) {
  PrsmReader reader(input_file_name);
  PrsmStrPtrVec prsm_str_ptrs;
  PrsmStrPtr prsm_str_ptr = reader.readOnePrsmStr();
  while (prsm_str_ptr != nullptr) {
    prsm_str_ptrs.push_back(prsm_str_ptr);
    prsm_str_ptr = reader.readOnePrsmStr();
  }
  reader.close();
  return prsm_str_ptrs;
}

PrsmPtrVec PrsmReader::readAllPrsms(const std::string &prsm_file_name,
                                    const std::string &db_file_name,
                                    const ModPtrVec  &fix_mod_list) {
  FastaIndexReaderPtr fasta_reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);
  PrsmReader reader(prsm_file_name);
  PrsmPtrVec prsm_ptrs;
  PrsmPtr prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  while (prsm_ptr != nullptr) {
    prsm_ptrs.push_back(prsm_ptr);
    prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  }
  reader.close();
  return prsm_ptrs;
}

}  // namespace prot
