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


#include "feature/msalign_writer.hpp"

namespace prot {

void MsalignWriter::write(std::ofstream &file, MatchEnvPtrVec &envs,
                              MsHeaderPtr header_ptr) {
  file << "BEGIN IONS" << std::endl;
  file << "ID=" << header_ptr->getId() << std::endl;
  file << "SCANS=" << header_ptr->getScansString() << std::endl;
  file << "RETENTION_TIME=" << header_ptr->getRetentionTime() << std::endl;
  if (header_ptr->getActivationPtr() != nullptr) {
    file << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
  }
  if (header_ptr->getMsLevel() > 1) {
    file << "MS_ONE_ID=" << header_ptr->getMsOneId() << std::endl;
    file << "MS_ONE_SCAN=" << header_ptr->getMsOneScan() << std::endl;
    file << "PRECURSOR_MZ=" << header_ptr->getPrecMonoMz() << std::endl;
    file << "PRECURSOR_CHARGE=" << header_ptr->getPrecCharge() << std::endl;
    file << "PRECURSOR_MASS=" <<  header_ptr->getPrecMonoMass() << std::endl;
    file << "PRECURSOR_INTENSITY=" << header_ptr->getPrecInte() << std::endl;
    if (header_ptr->getFeatureId() >= 0) {
      file << "FEATURE_ID=" << header_ptr->getFeatureId() << std::endl;
      file << "FEATURE_INTENSITY=" << header_ptr->getFeatureInte() << std::endl;
    }
  }

  for (size_t i = 0; i < envs.size(); i++) {
    MatchEnvPtr env = envs[i];
    EnvelopePtr theo_env = env->getTheoEnvPtr();
    RealEnvPtr real_env = env->getRealEnvPtr();
    file << real_env->getMonoMass();
    file << "\t" << theo_env->compIntensitySum();
    file << "\t" << theo_env->getCharge();
    file << std::endl;
  }
  file << "END IONS" << std::endl;
  file << std::endl;
}

void MsalignWriter::write(std::ofstream &file, DeconvMsPtr ms_ptr) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  file << "BEGIN IONS" << std::endl;
  file << "ID=" << header_ptr->getId() << std::endl;
  file << "SCANS=" << header_ptr->getScansString() << std::endl;
  file << "RETENTION_TIME=" << header_ptr->getRetentionTime() << std::endl;
  if (header_ptr->getActivationPtr() != nullptr) {
    file << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
  }
  if (header_ptr->getMsLevel() > 1) {
    file << "MS_ONE_ID=" << header_ptr->getMsOneId() << std::endl;
    file << "MS_ONE_SCAN=" << header_ptr->getMsOneScan() << std::endl;
    file << "PRECURSOR_MZ=" << header_ptr->getPrecMonoMz() << std::endl;
    file << "PRECURSOR_CHARGE=" << header_ptr->getPrecCharge() << std::endl;
    file << "PRECURSOR_MASS=" <<  header_ptr->getPrecMonoMass() << std::endl;
    file << "PRECURSOR_INTENSITY=" << header_ptr->getPrecInte() << std::endl;
    if (header_ptr->getFeatureId() >= 0) {
      file << "FEATURE_ID=" << header_ptr->getFeatureId() << std::endl;
      file << "FEATURE_INTENSITY=" << header_ptr->getFeatureInte() << std::endl;
    }
  }

  for (size_t i = 0; i < ms_ptr->size(); i++) {
    DeconvPeakPtr peak_ptr = ms_ptr->getPeakPtr(i);
    file << peak_ptr->getPosition();
    file << "\t" << peak_ptr->getIntensity();
    file << "\t" << peak_ptr->getCharge();
    file << std::endl;
  }
  file << "END IONS" << std::endl;
  file << std::endl;
}

void MsalignWriter::write(std::ofstream &file, DeconvMsPtr ms_ptr, int mslevel) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  file << "BEGIN IONS" << std::endl;
  file << "ID=" << header_ptr->getId() << std::endl;
  file << "SCANS=" << header_ptr->getScansString() << std::endl;
  file << "RETENTION_TIME=" << header_ptr->getRetentionTime() << std::endl;
  if (header_ptr->getActivationPtr() != nullptr) {
    file << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
  }
  if (mslevel > 1) {
    file << "MS_ONE_ID=" << header_ptr->getMsOneId() << std::endl;
    file << "MS_ONE_SCAN=" << header_ptr->getMsOneScan() << std::endl;
    file << "PRECURSOR_MZ=" << header_ptr->getPrecMonoMz() << std::endl;
    file << "PRECURSOR_CHARGE=" << header_ptr->getPrecCharge() << std::endl;
    file << "PRECURSOR_MASS=" <<  header_ptr->getPrecMonoMass() << std::endl;
    file << "PRECURSOR_INTENSITY=" << header_ptr->getPrecInte() << std::endl;
    if (header_ptr->getFeatureId() >= 0) {
      file << "FEATURE_ID=" << header_ptr->getFeatureId() << std::endl;
      file << "FEATURE_INTENSITY=" << header_ptr->getFeatureInte() << std::endl;
    }
  }

  for (size_t i = 0; i < ms_ptr->size(); i++) {
    DeconvPeakPtr peak_ptr = ms_ptr->getPeakPtr(i);
    file << peak_ptr->getPosition();
    file << "\t" << peak_ptr->getIntensity();
    file << "\t" << peak_ptr->getCharge();
    file << std::endl;
  }
  file << "END IONS" << std::endl;
  file << std::endl;
}


}  // namespace prot
