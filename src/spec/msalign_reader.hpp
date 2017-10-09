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


#ifndef PROT_SPEC_MSALIGN_READER_HPP_
#define PROT_SPEC_MSALIGN_READER_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <set>

#include "base/logger.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"

namespace prot {

class MsAlignReader {
 public:
  MsAlignReader(const std::string &file_name, int group_spec_num,
                ActivationPtr act_ptr, const std::set<std::string> skip_list):
      file_name_(file_name),
      group_spec_num_(group_spec_num),
      activation_ptr_(act_ptr),
      skip_list_(skip_list) {
        input_.open(file_name.c_str(), std::ios::in);
        if (!input_.is_open()) {
          LOG_ERROR("msalign file  " << file_name << " does not exist.");
          throw "msalign file does not exist.";
        }
      }

  std::vector<std::string> readOneSpectrum();

  DeconvMsPtr getNextMs();

  std::vector<SpectrumSetPtr> getNextSpectrumSet(SpParaPtr sp_para_ptr);

  void close();

 private:
  std::string file_name_;
  int group_spec_num_;
  std::ifstream input_;
  std::vector<std::string> spectrum_str_vec_;
  int current_ = 0;
  ActivationPtr activation_ptr_;
  DeconvMsPtr deconv_ms_ptr_ = nullptr;
  std::set<std::string> skip_list_;

  void readNext();
};

typedef std::shared_ptr<MsAlignReader> MsAlignReaderPtr;

}  // namespace prot
#endif
