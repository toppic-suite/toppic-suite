//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include "spec/msalign_writer.hpp"

namespace prot {

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
