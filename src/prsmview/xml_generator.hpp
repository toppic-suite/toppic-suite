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


#ifndef PROT_XML_GENERATOR_HPP_
#define PROT_XML_GENERATOR_HPP_

#include <map>
#include <string>
#include <vector>

#include "xercesc/util/PlatformUtils.hpp"

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/xml_writer.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsmview/prsm_view_mng.hpp"
#include "prsmview/anno_view.hpp"

namespace prot {

class XmlGenerator {
 public:
  XmlGenerator(PrsmParaPtr prsm_para_ptr,
               const std::string &exec_dir,
               const std::string &input_file_name,
               const std::string &fname_suffix);
  void process();

 private:
  void outputPrsms();
  void outputProteoforms();
  void outputProteins();
  void outputAllProteins();
  void outputFileList();
  void splitBySpeciesId();
  void splitByProtId();

  std::string input_file_ext_;
  PrsmViewMngPtr mng_ptr_;
  AnnoViewPtr anno_view_ptr_;
  FastaIndexReaderPtr fasta_reader_ptr_;
  std::vector<int> species_ids_;
  std::vector<int> prot_ids_;
  int writer_block_size_;
  std::vector<ExtendMsPtrVec> extend_ms_vec2d_;
  std::vector<DeconvMsPtrVec> deconv_ms_vec2d_;
  std::map<int, size_t> spec_id_extend_ms_map_;
};

typedef std::shared_ptr<XmlGenerator> XmlGeneratorPtr;

}  // namespace prot

#endif /* PROT_XML_GENERATOR_HPP_ */
