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

#ifndef TOPPIC_PRSM_PRSM_XML_WRITER_HPP_
#define TOPPIC_PRSM_PRSM_XML_WRITER_HPP_

#include <fstream>
#include <string>

#include "prsm/prsm.hpp"
#include "prsm/prsm_str.hpp"

namespace toppic {

class PrsmXmlWriter {
 public:
  explicit PrsmXmlWriter(const std::string& file_name);

  ~PrsmXmlWriter();

  void close();

  void write(const PrsmStrPtr& prsm_str_ptr);

  void writeVector(const PrsmStrPtrVec& prsm_str_ptr_vec);

  void write(const PrsmPtr& prsm_ptr);

  // Writes the PrSM as the merge of per-thread files (PrsmStrMerge) wrote
  // it, through PrsmReader::readOnePrsmLines: every line of the XML string
  // trimmed, empty lines left out; a stage that writes its PrSMs in order
  // itself gives the merge's file this way (toppic_test experiment 150).
  void writeTrimmedLines(const PrsmPtr& prsm_ptr);

  // The XML write(prsm_ptr) writes, as a string, and writing such a string:
  // what a thread that computes PrSMs for this file can do without the
  // file, which is written later in its own order. toXmlString does not
  // touch the file.
  std::string toXmlString(const PrsmPtr& prsm_ptr);

  void writeXmlString(const std::string& str);

  void writeVector(const PrsmPtrVec& prsm_ptrs);

  void writeVector2D(const PrsmPtrVec2D& prsm_ptrs);

  void writeVector3D(const PrsmPtrVec3D& prsm_ptrs);

 private:
  std::ofstream file_;

  std::string file_name_;
};

using PrsmXmlWriterPtr = std::shared_ptr<PrsmXmlWriter>;
using PrsmXmlWriterPtrVec = std::vector<PrsmXmlWriterPtr>;

} /* namespace toppic */

#endif /* PRSM_WRITER_HPP_ */
