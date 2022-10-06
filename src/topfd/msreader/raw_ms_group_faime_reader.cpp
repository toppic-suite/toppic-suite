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
#include "ms/env/real_env.hpp"
#include "ms/env/prec_env.hpp"
#include "topfd/msreader/raw_ms_group_faime_reader.hpp" 

namespace toppic {

RawMsGroupFaimeReader::RawMsGroupFaimeReader(const std::string & file_name, 
                                             bool missing_level_one,
                                             std::string activation,
                                             double isolation_window,
                                             int fraction_id) {

  reader_ptr_ = std::make_shared<PwMsReader>(file_name, activation, 
                                             isolation_window);
  missing_level_one_ = missing_level_one;
  fraction_id_ = fraction_id;
  if (!missing_level_one_) {
    RawMsPtr ms_ptr = readNextRawMs(); 
    if (ms_ptr == nullptr) {
      LOG_ERROR("The file " << file_name << " does not contain spectra!");
      exit(EXIT_FAILURE);
    }
    int idx = 0;
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
      ms_ptr = readNextRawMs();
    }
    reader_ptr_ = std::make_shared<PwMsReader>(file_name, activation, 
                                               isolation_window);
    total_ms_one_num_ = idx;
    cur_ms_one_idx_ = 0;
  }
}

RawMsPtr RawMsGroupFaimeReader::readNextRawMs() {
  reader_ptr_->readNext();
  PeakPtrVec peak_list = reader_ptr_->getPeakList();
  MsHeaderPtr header_ptr = reader_ptr_->getHeaderPtr();
  if (header_ptr == nullptr) {
    return nullptr;
  }
  header_ptr->setFractionId(fraction_id_);
  RawMsPtr ms_ptr = std::make_shared<Ms<PeakPtr> >(header_ptr, peak_list);
  return ms_ptr;
}

RawMsGroupPtr RawMsGroupFaimeReader::getNextMsGroupPtrWithFaime() {
  if (missing_level_one_) {
    RawMsPtr null_ms_one_ptr(nullptr);
    RawMsPtr ms_two_ptr = readNextRawMs();
    if (ms_two_ptr == nullptr) {
      return nullptr;
    }
    RawMsPtrVec ms_two_ptr_vec;
    ms_two_ptr_vec.push_back(ms_two_ptr);
    RawMsGroupPtr ms_group_ptr = std::make_shared<RawMsGroup>(null_ms_one_ptr, ms_two_ptr_vec);
    return ms_group_ptr;
  }
  while (true) {
    while (cur_ms_one_idx_ < total_ms_one_num_) {
      int ms1_scan = ms_one_scans_[cur_ms_one_idx_];
      // if the ms1 scan to be deconvolted has not been read
      if (cur_last_ms_two_scan_map_.count(ms1_scan) == 0) {
        break;
      }
      if (cur_last_ms_two_scan_map_[ms1_scan] < last_ms_two_scan_map_[ms1_scan]) {
        break;
      }
      // if all ms1 scan and all ms2 scan has been read
      if (cur_last_ms_two_scan_map_[ms1_scan] >= last_ms_two_scan_map_[ms1_scan]) {
        if (last_ms_two_scan_map_[ms1_scan] == -1) {
          // if the ms1 scan does not have ms2, return a ms group with the ms1
          // scan only 
          RawMsGroupPtr ms_group_ptr = std::make_shared<RawMsGroup>(ms_one_ptr_map_[ms1_scan], RawMsPtrVec());
          ms_one_ptr_map_[ms1_scan] = nullptr;
          cur_ms_one_idx_ = cur_ms_one_idx_ + 1;
          return ms_group_ptr;
        }
        else {
          RawMsGroupPtr ms_group_ptr = std::make_shared<RawMsGroup>(ms_one_ptr_map_[ms1_scan], ms_two_ptr_vec_map_[ms1_scan]);
          ms_one_ptr_map_[ms1_scan] = nullptr;
          ms_two_ptr_vec_map_[ms1_scan].clear();
          cur_ms_one_idx_ = cur_ms_one_idx_ + 1;
          return ms_group_ptr;
        }
      }
    }

    RawMsPtr ms_ptr = readNextRawMs();
    if (ms_ptr == nullptr) {
      break;
    }
    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    if (header_ptr->getMsLevel() == 1) {
      int ms1_scan = header_ptr->getFirstScanNum();
      ms_one_ptr_map_[ms1_scan] = ms_ptr;
      cur_last_ms_two_scan_map_[ms1_scan] = -1;
    }
    else {
      int ms2_scan = header_ptr->getFirstScanNum();
      int ms1_scan = header_ptr->getMsOneScan();
      if (ms_two_ptr_vec_map_.count(ms1_scan) == 0) {
        ms_two_ptr_vec_map_[ms1_scan] = RawMsPtrVec();
      }
      ms_two_ptr_vec_map_[ms1_scan].push_back(ms_ptr);
      cur_last_ms_two_scan_map_[ms1_scan] = ms2_scan;
    }
  }
  return nullptr;
}

// refine precursor charge and mz 
MatchEnvPtr refinePrecChrgFaime(RawMsPtr ms_one, RawMsPtr ms_two, 
                                double max_mass, int max_charge) {
  MsHeaderPtr header_two = ms_two->getMsHeaderPtr();
  double prec_win_begin = header_two->getPrecWinBegin();
  double prec_win_end = header_two->getPrecWinEnd();

  PeakPtrVec peak_list = ms_one->getPeakPtrVec();
  LOG_DEBUG("start refine precursor " << " peak num " << peak_list.size());
  MatchEnvPtr match_env_ptr = prec_env::deconv(prec_win_begin, prec_win_end, peak_list,  
                                               max_mass, max_charge);
  if (match_env_ptr != nullptr) {
    RealEnvPtr env_ptr = match_env_ptr->getRealEnvPtr();
    header_two->setPrecMonoMz(env_ptr->getMonoMz());
    header_two->setPrecCharge(env_ptr->getCharge());
    header_two->setPrecInte(env_ptr->compIntensitySum());
    LOG_DEBUG("prec mz " << env_ptr->getMonoMz() << " prec charge " << env_ptr->getCharge());
  } else {
    header_two->setPrecMonoMz(0);
    header_two->setPrecCharge(0);
    header_two->setPrecInte(0);
    LOG_INFO("EMPTY ENVELOPE POINTER");
  }
  return match_env_ptr;
}

void RawMsGroupFaimeReader::obtainPrecEnvs(RawMsGroupPtr ms_group_ptr, 
                                           MatchEnvPtrVec &env_ptr_vec,
                                           double max_mass, int max_charge) {
  RawMsPtr ms_one_ptr = ms_group_ptr->getMsOnePtr();
  RawMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();

  for (size_t i = 0; i < ms_two_ptr_vec.size(); i++) {
    RawMsPtr ms_two_ptr = ms_two_ptr_vec[i];
    MatchEnvPtr match_env_ptr = refinePrecChrgFaime(ms_one_ptr, ms_two_ptr, 
                                                    max_mass, max_charge);
    if (match_env_ptr != nullptr) {
      env_ptr_vec.push_back(match_env_ptr);
    }
  }
}


}
