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


#ifndef PROT_LOCAL_MNG_HPP_
#define PROT_LOCAL_MNG_HPP_


#include <string>
#include <vector>
#include <utility>

#include "prsm/prsm_para.hpp"
#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"

namespace prot {

const int LEFT_SUP_LIMIT = 10;
const int RIGHT_SUP_LIMIT = 10;
const int DESC_MATCH_LIMIT = 5;

typedef std::pair<PtmPtr, PtmPtr> PtmPair;
typedef std::vector<PtmPair> PtmPairVec;

class LocalMng {
 public:
  LocalMng(PrsmParaPtr prsm_para_ptr, const std::string& local_threshold,
           const std::string& residueModFileName, double max_ptm_mass,
           const std::string &input_file_ext, const std::string &output_file_ext):
      prsm_para_ptr_(prsm_para_ptr),
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext),
      residueModFileName_(residueModFileName),
      threshold_(std::stod(local_threshold)),
      theta_(0.994),
      beta_(0.8),
      min_mass_(prsm_para_ptr->getSpParaPtr()->getMinMass()),
      max_ptm_mass_(max_ptm_mass),
      p1_(0.915258),
      p2_(21.1822) {}

  PrsmParaPtr prsm_para_ptr_;

  std::string input_file_ext_;
  std::string output_file_ext_;
  std::string residueModFileName_;

  double threshold_, theta_, beta_;
  double min_mass_, max_ptm_mass_;
  double p1_, p2_;
};

typedef std::shared_ptr<LocalMng> LocalMngPtr;

}  // namespace prot

#endif /* PROT_LOCAL_MNG_HPP_ */
