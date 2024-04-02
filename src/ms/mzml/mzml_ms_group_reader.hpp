//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_MS_MZML_MZML_MS_GROUP_FAIME_READER_HPP_
#define TOPPIC_MS_MZML_MZML_MS_GROUP_FAIME_READER_HPP_

#include "ms/env/match_env.hpp"
#include "ms/mzml/mzml_ms_group.hpp"
#include "ms/mzml/pw_ms_reader.hpp"

namespace toppic {

class MzmlMsGroupReader {
 public:
  MzmlMsGroupReader(const std::string & file_name, 
                    double isolation_window, 
                    std::string activation, 
                    int fraction_id,
                    bool is_faims,
                    double faims_voltage,
                    bool missing_level_one); 

  MzmlMsGroupPtr getNextMsGroupPtr();

  PwMsReaderPtr getReaderPtr() {return reader_ptr_;}
  int getInputSpNum() {return reader_ptr_->getInputSpNum();}
  bool checkCentroidData() {return reader_ptr_->checkCentroidData();}

  void getMs1Map(PeakPtrVec2D &ms1_mzml_peaks, MsHeaderPtr2D &ms2_header_ptr_2d);
  void getMs2Map(PeakPtrVec2D &ms2_mzml_peaks, double base_mz);

 private:
  PwMsReaderPtr reader_ptr_;
  int fraction_id_ = 0;
  bool is_faims_ = false;
  double faims_voltage_ = -1;
  bool missing_level_one_ = false;

  std::map<int, MzmlMsPtr> ms_one_ptr_map_;
  std::map<int, int> ms_one_scan_idx_map_;
  std::map<int, MzmlMsPtrVec> ms_two_ptr_vec_map_;
  std::map<int, int> last_ms_two_scan_map_;
  std::vector<int> ms_one_scans_;

  int total_ms_one_num_;
  int cur_ms_one_idx_;
  std::map<int, int> cur_last_ms_two_scan_map_;

  int ms_one_cnt_ = 0;
  int ms_two_cnt_ = 0;

  MzmlMsPtr readNextMzmlMs();

  void initMs2Ms1Map();

  MzmlMsPtr readNextMs2MzmlMs(); 

  MzmlMsGroupPtr getMs2OnlyMsGroupPtr();

  MzmlMsGroupPtr getMs1Ms2MsGroupPtr(); 
};

typedef std::shared_ptr<MzmlMsGroupReader> MzmlMsGroupReaderPtr;

}

#endif
