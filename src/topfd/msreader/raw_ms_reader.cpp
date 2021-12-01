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
#include "topfd/msreader/raw_ms_reader.hpp" 

namespace toppic {

RawMsReader::RawMsReader(const std::string & file_name, double isolation_window) {
  reader_ptr_ = std::make_shared<PwMsReader>(file_name, isolation_window);
}

RawMsReader::RawMsReader(const std::string & file_name, const std::string & activation,
                         double isolation_window) {
  reader_ptr_ = std::make_shared<PwMsReader>(file_name, activation, isolation_window);
}

RawMsPtr RawMsReader::getNextMs(double max_mass, int max_charge) {
  reader_ptr_->readNext();
  PeakPtrVec peak_list = reader_ptr_->getPeakList();
  MsHeaderPtr header_ptr = reader_ptr_->getHeaderPtr();
  if (header_ptr == nullptr) {
    return nullptr;
  }
  RawMsPtr ms_ptr = std::make_shared<Ms<PeakPtr> >(header_ptr, peak_list);

  // update ms_1 
  if (header_ptr->getMsLevel() == 1) {
    ms_one_ = ms_ptr;
  }
  if (header_ptr->getMsLevel() == 2 && ms_one_ != nullptr) {
    header_ptr->setMsOneId(ms_one_->getMsHeaderPtr()->getId());
    header_ptr->setMsOneScan(ms_one_->getMsHeaderPtr()->getFirstScanNum());
  }
  if (do_refine_prec_mass_ && header_ptr->getMsLevel() == 2 && ms_one_ != nullptr) {
    refinePrecChrg(ms_one_, ms_ptr, max_mass, max_charge);
  } else {
    if (header_ptr->getPrecSpMz() != 0.0) {
      header_ptr->setPrecMonoMz(header_ptr->getPrecSpMz());
    }
  }
  return ms_ptr;
}

void RawMsReader::getMs1Peaks(PeakPtrVec2D &raw_peaks, double cur_voltage) {
  while (true) {
    reader_ptr_->readNext();
    MsHeaderPtr header_ptr = reader_ptr_->getHeaderPtr();
    if (header_ptr == nullptr) {
      break;
    }
    if (header_ptr->getMsLevel() == 1 && header_ptr->getVoltage() == -1) {//if not a FAIME data
      PeakPtrVec peak_list = reader_ptr_->getPeakList();
      raw_peaks.push_back(peak_list);
    }
    else if (header_ptr->getMsLevel() == 1 && header_ptr->getVoltage() == cur_voltage) {//FAIME data
      PeakPtrVec peak_list = reader_ptr_->getPeakList();
      raw_peaks.push_back(peak_list);
    }
  }
}

// refine precursor charge and mz 
void RawMsReader::refinePrecChrg(RawMsPtr ms_one, RawMsPtr ms_two, 
                                 double max_mass, int max_charge) {
  MsHeaderPtr header_two = ms_two->getMsHeaderPtr();
  double prec_win_begin = header_two->getPrecWinBegin();
  double prec_win_end = header_two->getPrecWinEnd();
  //int prec_charge = header_two->getPrecCharge();

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
    LOG_DEBUG("EMPTY ENVELOPE POINTER");
  }
}

}
