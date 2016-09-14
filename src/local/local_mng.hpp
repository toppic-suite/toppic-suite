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


#ifndef PROT_LOCAL_MNG_HPP_
#define PROT_LOCAL_MNG_HPP_

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
           const std::string &input_file_ext, const std::string &output_file_ext);

  PrsmParaPtr prsm_para_ptr_;

  std::string residueModFileName_;
  std::string input_file_ext_;
  std::string output_file_ext_;

  double threshold_, theta_, beta_;
  double min_mass_, max_ptm_mass_;
  double p1_, p2_;
};

typedef std::shared_ptr<LocalMng> LocalMngPtr;

inline LocalMng::LocalMng(PrsmParaPtr prsm_para_ptr, const std::string& local_threshold, 
                          const std::string& residueModFileName, double max_ptm_mass,
                          const std::string &input_file_ext, const std::string &output_file_ext) {

  prsm_para_ptr_ = prsm_para_ptr;
  input_file_ext_ = input_file_ext;
  output_file_ext_ = output_file_ext;
  residueModFileName_ = residueModFileName;
  min_mass_ = prsm_para_ptr->getSpParaPtr()->getMinMass();
  threshold_ = std::stod(local_threshold);
  max_ptm_mass_ = max_ptm_mass;
  theta_ = 0.994;
  beta_ = 0.8;
  p1_ = 0.915258;
  p2_ = 21.1822;
}

}

#endif /* PROT_LOCAL_MNG_HPP_ */
