//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_MS_MZML_PW_MS_READER_HPP_
#define TOPPIC_MS_MZML_PW_MS_READER_HPP_

#include <memory>
#include <string>

#include "pwiz/data/msdata/DefaultReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"

#include "ms/mzml/mzml_ms.hpp"
#include "ms/mzml/mzml_profile.hpp"

namespace toppic {

typedef std::shared_ptr<pwiz::msdata::MSDataFile> MSDataFilePtr;

class PwMsReader {
 public:
  explicit PwMsReader(const std::string & file_name); 

  explicit PwMsReader(const std::string & file_name, double isolation_window);
  
  explicit PwMsReader(const std::string & file_name, double isolation_window,
                      std::string activation);
  int readNext();
  int readNextWithVoltage(double voltage); 
  MsHeaderPtr getHeaderPtr() {return header_ptr_;}
  PeakPtrVec getPeakList() {return peak_list_;}
  int getInputSpNum() {return input_sp_num_;}
  bool checkCentroidData(); 

  MzmlProfilePtr readProfile();

  // reset indexes
  void resetIndexes(); 

 private:
  std::string file_name_;
  std::string activation_;
  double isolation_window_;
  //functions for waters instrument is different from others 
  bool is_waters_instrument_ = false;

  int input_sp_num_;
  int input_sp_id_ = 0;
  int ms1_cnt_ = 0;
  int ms2_cnt_ = 0;
  int prev_ms1_scan_id_ = -1;
  PeakPtrVec peak_list_;
  MsHeaderPtr header_ptr_;

  // pwiz reader
  pwiz::msdata::DefaultReaderList readers_;
  MSDataFilePtr msd_ptr_;
  pwiz::msdata::SpectrumListPtr spec_list_ptr_;

  void init(const std::string &file_name);

  bool readOneMs(int sp_id, PeakPtrVec &peak_list, MsHeaderPtr &header_ptr);

  bool checkWatersInstrument(); 

  PeakPtrVec parsePeaks(pwiz::msdata::SpectrumPtr cur_spec_ptr); 

  int parseNum(std::string &id, int default_scan); 

  void parseScanNum(MsHeaderPtr header_ptr, 
                    pwiz::msdata::SpectrumInfo &spec_info); 

  void parsePrecursor(MsHeaderPtr header_ptr, 
                      pwiz::msdata::SpectrumInfo &spec_info,
                      pwiz::msdata::SpectrumPtr cur_spec_ptr);

  void parseActivation(MsHeaderPtr header_ptr, 
                       pwiz::msdata::SpectrumInfo &spec_info,
                       pwiz::msdata::SpectrumPtr cur_spec_ptr); 

  double parseFaims(pwiz::msdata::SpectrumPtr cur_spec_ptr);
};

typedef std::shared_ptr<PwMsReader> PwMsReaderPtr;

}  // namespace toppic

#endif
