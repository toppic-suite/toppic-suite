//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

  std::vector<int> cluster_ids_;

  std::vector<int> prot_ids_;

  int writer_block_size_;

  std::vector<ExtendMsPtrVec> extend_ms_vec2d_;

  std::vector<DeconvMsPtrVec> deconv_ms_vec2d_;

  std::map<int, size_t> spec_id_extend_ms_map_;
};

typedef std::shared_ptr<XmlGenerator> XmlGeneratorPtr;

}  // namespace prot

#endif /* PROT_XML_GENERATOR_HPP_ */
