//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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


#ifndef TOPPIC_TOPFD_MSREADER_PW_MS_READER_HPP_
#define TOPPIC_TOPFD_MSREADER_PW_MS_READER_HPP_

#include <memory>
#include <string>

#include "pwiz/data/msdata/DefaultReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"

#include "ms/spec/raw_ms.hpp"

namespace toppic {

typedef std::shared_ptr<pwiz::msdata::MSDataFile> MSDataFilePtr;

class PwMsReader {
 public:
  explicit PwMsReader(const std::string & file_name);
  int readNext();
  PeakPtrVec getPeakList() {return peak_list_;}
  MsHeaderPtr getHeaderPtr() {return header_ptr_;}
  int getInputSpNum() {return input_sp_num_;}
  bool readOneMs(int sp_id, PeakPtrVec &peak_list, MsHeaderPtr &header_ptr);

 private:
  std::string file_name_;
  int input_sp_num_;
  int input_sp_id_;
  int output_sp_id_;

  int ms1_cnt = 0;
  int ms2_cnt = 0;
  PeakPtrVec peak_list_;
  MsHeaderPtr header_ptr_;

  // pwiz reader
  pwiz::msdata::DefaultReaderList readers_;
  MSDataFilePtr msd_ptr_;
  pwiz::msdata::SpectrumListPtr spec_list_ptr_;
};

typedef std::shared_ptr<PwMsReader> PwMsReaderPtr;

}  // namespace toppic

#endif
