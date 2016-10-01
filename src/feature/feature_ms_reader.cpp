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


#include "base/logger.hpp"
#include "feature/real_env.hpp"
#include "feature/prec_env.hpp"
#include "feature/feature_ms_reader.hpp" 

namespace prot {

FeatureMsReader::FeatureMsReader(std::string &file_name) {
  reader_ptr_ = RawMsReaderPtr(new RawMsReader(file_name));
}

RawMsPtr FeatureMsReader::getNextMs(double prec_win_size) {
  reader_ptr_->readNext();
  PeakPtrVec peak_list = reader_ptr_->getPeakList();
  MsHeaderPtr header_ptr = reader_ptr_->getHeaderPtr();
  if (header_ptr == nullptr) {
    return nullptr;
  }
  RawMsPtr ms_ptr(new Ms<PeakPtr>(header_ptr, peak_list));

  // update ms_1 
  if (header_ptr->getMsLevel() == 1) {
    ms_one_ = ms_ptr;
  }
  if (header_ptr->getMsLevel() == 2 && ms_one_ != nullptr) {
    header_ptr->setMsOneId(ms_one_->getMsHeaderPtr()->getId());
    header_ptr->setMsOneScan(ms_one_->getMsHeaderPtr()->getFirstScanNum());
  }
  if (do_refine_prec_mass_ && header_ptr->getMsLevel() == 2 && ms_one_ != nullptr) {
    refinePrecChrg(ms_one_, ms_ptr, prec_win_size);
  } else {
    if (header_ptr->getPrecSpMz() != 0.0) {
      header_ptr->setPrecMonoMz(header_ptr->getPrecSpMz());
    }
  }
  return ms_ptr;
}

// refine precursor charge and mz 
void FeatureMsReader::refinePrecChrg(RawMsPtr ms_one, RawMsPtr ms_two, 
                                     double prec_win_size) {
  MsHeaderPtr header_two = ms_two->getMsHeaderPtr();
  double prec_avg_mz = header_two->getPrecSpMz();
  int prec_charge = header_two->getPrecCharge();

  PeakPtrVec peak_list = ms_one->getPeakPtrVec();
  LOG_DEBUG("start refine precursor " << " peak num " << peak_list.size());
  RealEnvPtr env_ptr = PrecEnv::deconv(prec_win_size, peak_list, prec_avg_mz, 
                                       prec_charge);
  if (env_ptr != nullptr) {
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
