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


#include "common/util/logger.hpp"
#include "ms/env/exp_env.hpp"
#include "ms/mzml/mzml_ms_group_reader.hpp" 

#include "topfd/deconv/prec_env.hpp"

namespace toppic {

MzmlMsGroupReader::MzmlMsGroupReader(const std::string & file_name, 
                                     double isolation_window, 
                                     std::string activation, 
                                     int fraction_id,
                                     double faims_voltage,
                                     bool missing_level_one) {
  fraction_id_ = fraction_id_;
  faims_voltage_ = faims_voltage;
  missing_level_one_ = missing_level_one;
  reader_ptr_ = std::make_shared<PwMsReader>(file_name, 
                                             isolation_window, 
                                             activation);
  MzmlMsPtr ms_ptr = readNextMzmlMs(); 
  if (ms_ptr == nullptr) {
    LOG_ERROR("The file " << file_name << " does not contain spectra!");
    exit(EXIT_FAILURE);
  }

  if (!missing_level_one_) {
    reader_ptr_->resetIndexes();
    initMs2Ms1Map();
    reader_ptr_->resetIndexes();
    cur_ms_one_idx_ = 0;
  }
}

MzmlMsPtr MzmlMsGroupReader::readNextMzmlMs() {
  //Read ms with a specific voltage
  if (faims_voltage_ < 0) {
    reader_ptr_->readNext();
  }
  else {
    reader_ptr_->readNextWithVoltage(faims_voltage_);
  }
  MsHeaderPtr header_ptr = reader_ptr_->getHeaderPtr();
  if (header_ptr == nullptr) {
    return nullptr;
  }
  header_ptr->setFractionId(fraction_id_);
  PeakPtrVec peak_list = reader_ptr_->getPeakList();
  MzmlMsPtr ms_ptr = std::make_shared<Ms<PeakPtr> >(header_ptr, peak_list);
  return ms_ptr;
}

MzmlMsPtr MzmlMsGroupReader::readNextMs2MzmlMs() {
  while (true) {
    MzmlMsPtr ms_ptr = readNextMzmlMs();
    if (ms_ptr == nullptr) {
      return nullptr;
    }
    if (ms_ptr->getMsHeaderPtr()->getMsLevel() == 2) {
      return ms_ptr;
    }
  }
}

void MzmlMsGroupReader::initMs2Ms1Map() {
  int idx = 0;
  MzmlMsPtr ms_ptr = readNextMzmlMs(); 
  while (ms_ptr != nullptr) {
    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    if (header_ptr->getMsLevel() == 1) {
      int ms1_scan = header_ptr->getFirstScanNum();
      //set inital value of the last ms two scan of the ms1 scan to -1 
      last_ms_two_scan_map_[ms1_scan] = -1;
      ms_one_scans_.push_back(ms1_scan);
      ms_one_scan_idx_map_[ms1_scan] = idx;
      idx = idx + 1;
    }
    else if (header_ptr->getMsLevel() == 2) {
      int ms2_scan = header_ptr->getFirstScanNum();
      int ms1_scan = header_ptr->getMsOneScan();
      last_ms_two_scan_map_[ms1_scan] = ms2_scan;
    }
    ms_ptr = readNextMzmlMs();
  }
  total_ms_one_num_ = idx;
}


MzmlMsGroupPtr MzmlMsGroupReader::getMs2OnlyMsGroupPtr() {
  MzmlMsPtr null_ms_one_ptr(nullptr);
  MzmlMsPtr ms_two_ptr = readNextMs2MzmlMs();
  if (ms_two_ptr == nullptr) {
    return nullptr;
  }
  MzmlMsPtrVec ms_two_ptr_vec;
  ms_two_ptr_vec.push_back(ms_two_ptr);
  MzmlMsGroupPtr ms_group_ptr = std::make_shared<MzmlMsGroup>(null_ms_one_ptr, ms_two_ptr_vec);
  return ms_group_ptr;
}

MzmlMsGroupPtr MzmlMsGroupReader::getMs1Ms2MsGroupPtr() {
  if (cur_ms_one_idx_ >= total_ms_one_num_) {
    return nullptr;
  }
  int ms1_scan = ms_one_scans_[cur_ms_one_idx_];
  // find the ms group for ms1_scan
  while (true) {
    // if all ms1 scan and all ms2 scan has been read, then return new ms group
    if (cur_last_ms_two_scan_map_.count(ms1_scan) > 0 &&
        cur_last_ms_two_scan_map_[ms1_scan] >= last_ms_two_scan_map_[ms1_scan]) {
      if (last_ms_two_scan_map_[ms1_scan] == -1) {
        // if the ms1 scan does not have ms2, return a ms group with the ms1
        // scan only 
        MzmlMsGroupPtr ms_group_ptr = std::make_shared<MzmlMsGroup>(ms_one_ptr_map_[ms1_scan], MzmlMsPtrVec());
        ms_one_ptr_map_[ms1_scan] = nullptr;
        cur_ms_one_idx_ = cur_ms_one_idx_ + 1;
        return ms_group_ptr;
      }
      else {
        MzmlMsPtrVec ms2_vec = ms_two_ptr_vec_map_[ms1_scan];
        // add ms one spec id for ms two later in deconv_process
        MzmlMsGroupPtr ms_group_ptr 
          = std::make_shared<MzmlMsGroup>(ms_one_ptr_map_[ms1_scan], ms_two_ptr_vec_map_[ms1_scan]);
        ms_one_ptr_map_[ms1_scan] = nullptr;
        ms_two_ptr_vec_map_[ms1_scan].clear();
        cur_ms_one_idx_ = cur_ms_one_idx_ + 1;
        return ms_group_ptr;
      }
    }
    // otherwise read more scans
    else {
      MzmlMsPtr ms_ptr = readNextMzmlMs();
      if (ms_ptr == nullptr) {
        break;
      }
      MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
      if (header_ptr->getMsLevel() == 1) {
        int cur_ms1_scan = header_ptr->getFirstScanNum();
        ms_one_ptr_map_[cur_ms1_scan] = ms_ptr;
        cur_last_ms_two_scan_map_[cur_ms1_scan] = -1;
      }
      else {
        int cur_ms1_scan = header_ptr->getMsOneScan();
        int cur_ms2_scan = header_ptr->getFirstScanNum();
        if (ms_two_ptr_vec_map_.count(cur_ms1_scan) == 0) {
          ms_two_ptr_vec_map_[cur_ms1_scan] = MzmlMsPtrVec();
        }
        ms_two_ptr_vec_map_[cur_ms1_scan].push_back(ms_ptr);
        cur_last_ms_two_scan_map_[cur_ms1_scan] = cur_ms2_scan;
      }
    }
  }
  // this line should not be reached
  return nullptr;
}

MzmlMsGroupPtr MzmlMsGroupReader::getNextMsGroupPtr() {
  if (missing_level_one_) {
    return getMs2OnlyMsGroupPtr();
  }
  else {
    return getMs1Ms2MsGroupPtr();
  }
}

}
