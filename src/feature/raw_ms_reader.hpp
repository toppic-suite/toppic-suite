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


#ifndef PROT_FEATURE_RAW_MS_READER_HPP_
#define PROT_FEATURE_RAW_MS_READER_HPP_

#include <memory>
#include <vector>
#include <string>

#include "pwiz/data/msdata/DefaultReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"

#include "spec/raw_ms.hpp"

namespace prot {

typedef std::shared_ptr<pwiz::msdata::MSDataFile> MSDataFilePtr;

class RawMsReader {
 public:
  RawMsReader(std::string &file_name);
  int readNext();
  PeakPtrVec getPeakList() {return peak_list_;}
  MsHeaderPtr getHeaderPtr() {return header_ptr_;}
  int getInputSpNum() {return input_sp_num_;}

 private:
  std::string file_name_;
  int input_sp_num_;
  int input_sp_id_;
  int output_sp_id_;

  int ms1_cnt = 0;
  int ms2_cnt = 0;
  PeakPtrVec peak_list_;
  MsHeaderPtr header_ptr_;
  
  //pwiz reader
  pwiz::msdata::DefaultReaderList readers_;
	MSDataFilePtr msd_ptr_; 
  pwiz::msdata::SpectrumListPtr spec_list_ptr_;
};

typedef std::shared_ptr<RawMsReader> RawMsReaderPtr;

}

#endif
