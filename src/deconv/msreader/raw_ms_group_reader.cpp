//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
#include "deconv/env/real_env.hpp"
#include "deconv/env/prec_env.hpp"
#include "deconv/msreader/raw_ms_group_reader.hpp" 

namespace toppic {

RawMsGroupReader::RawMsGroupReader(const std::string & file_name, bool missing_level_one) {
  reader_ptr_ = std::make_shared<PwMsReader>(file_name);
  missing_level_one_ = missing_level_one;
  if (!missing_level_one_) {
    ms_one_ptr_ = readNextRawMs();
    if (ms_one_ptr_ == nullptr) {
      LOG_ERROR("The file " << file_name << " does not contain spectra!");
      exit(EXIT_FAILURE);
    }
    if (ms_one_ptr_->getMsHeaderPtr()->getMsLevel() != 1) {
      LOG_ERROR("The first spectrum in " << file_name << " is not an MS1 spectrum!");
      exit(EXIT_FAILURE);
    }
  }
}

RawMsPtr RawMsGroupReader::readNextRawMs() {
  reader_ptr_->readNext();
  PeakPtrVec peak_list = reader_ptr_->getPeakList();
  MsHeaderPtr header_ptr = reader_ptr_->getHeaderPtr();
  if (header_ptr == nullptr) {
    return nullptr;
  }
  RawMsPtr ms_ptr = std::make_shared<Ms<PeakPtr> >(header_ptr, peak_list);
  return ms_ptr;
}

RawMsGroupPtr RawMsGroupReader::getNextMsGroupPtr() {
  if (missing_level_one_) {
    RawMsPtr null_ms_one_ptr(nullptr);
    RawMsPtrVec ms_two_ptr_vec;
    RawMsPtr ms_two_ptr = readNextRawMs();
    if (ms_two_ptr == nullptr) {
      return nullptr;
    }
    ms_two_ptr_vec.push_back(ms_two_ptr);
    RawMsGroupPtr ms_group_ptr = std::make_shared<RawMsGroup>(null_ms_one_ptr, ms_two_ptr_vec);
    return ms_group_ptr;
  }

  // if not missing level one
  if (ms_one_ptr_ == nullptr) {
    return nullptr;
  }

  RawMsPtrVec ms_two_ptr_vec;
  RawMsPtr new_ms_one_ptr;
  RawMsPtr ms_ptr;
  while ((ms_ptr = readNextRawMs()) != nullptr) {
    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    if (header_ptr->getMsLevel() == 1) {
      new_ms_one_ptr = ms_ptr;
      break;
    }
    else {
      header_ptr->setMsOneId(ms_one_ptr_->getMsHeaderPtr()->getId());
      header_ptr->setMsOneScan(ms_one_ptr_->getMsHeaderPtr()->getFirstScanNum());
      ms_two_ptr_vec.push_back(ms_ptr);
    }
  }
  RawMsGroupPtr ms_group_ptr = std::make_shared<RawMsGroup>(ms_one_ptr_, ms_two_ptr_vec);
  ms_one_ptr_ = new_ms_one_ptr;
  return ms_group_ptr;
}

// refine precursor charge and mz 
MatchEnvPtr refinePrecChrg(RawMsPtr ms_one, RawMsPtr ms_two, 
                           double prec_win_size, int max_charge) {
  MsHeaderPtr header_two = ms_two->getMsHeaderPtr();
  double prec_avg_mz = header_two->getPrecSpMz();
  int prec_charge = header_two->getPrecCharge();

  PeakPtrVec peak_list = ms_one->getPeakPtrVec();
  LOG_DEBUG("start refine precursor " << " peak num " << peak_list.size());
  MatchEnvPtr match_env_ptr = prec_env::deconv(prec_win_size, peak_list, prec_avg_mz, 
                                               prec_charge, max_charge);
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
    LOG_DEBUG("EMPTY ENVELOPE POINTER");
  }
  return match_env_ptr;
}

void RawMsGroupReader::obtainPrecEnvs(RawMsGroupPtr ms_group_ptr, 
                                      MatchEnvPtrVec &env_ptr_vec,
                                      double prec_win_size, int max_charge) {
  RawMsPtr ms_one_ptr = ms_group_ptr->getMsOnePtr();
  RawMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();
  for (size_t i = 0; i < ms_two_ptr_vec.size(); i++) {
    RawMsPtr ms_two_ptr = ms_two_ptr_vec[i];
    MatchEnvPtr match_env_ptr = refinePrecChrg(ms_one_ptr, ms_two_ptr, prec_win_size, max_charge);
    if (match_env_ptr != nullptr) {
      env_ptr_vec.push_back(match_env_ptr);
      //remove precursor peaks
      RealEnvPtr env_ptr = match_env_ptr->getRealEnvPtr();
      PeakPtrVec peaks = ms_one_ptr->getPeakPtrVec();
      for (int p = 0; p < env_ptr->getPeakNum(); p++) {
        if (env_ptr->isExist(p)) {
          int idx = env_ptr->getPeakIdx(p);
          peaks[idx]->setIntensity(0);
        }
      }
    }
  }
}


}
