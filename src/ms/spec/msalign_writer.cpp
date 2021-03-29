//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this output_ except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#include <iomanip>

#include "common/util/logger.hpp"
#include "ms/spec/msalign_writer.hpp"

namespace toppic {

MsAlignWriter::MsAlignWriter(const std::string &file_name) {
      output_.open(file_name);
      if (!output_.is_open()) {
        LOG_ERROR("Can not open the msalign file  " << file_name << "!");
        exit(EXIT_FAILURE);
      }
      output_.precision(16);
    }

MsAlignWriter::~MsAlignWriter() {
  if (output_.is_open()) {
    output_.close();
  }
}

void MsAlignWriter::close() {
  if (output_.is_open()) {
    output_.close();
  }
}

void MsAlignWriter::writePara(const std::string &para_str) {
  output_ << para_str << "\n";
}

void MsAlignWriter::write(DeconvMsPtr ms_ptr) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  output_ << std::fixed;
  output_ << "BEGIN IONS" << std::endl;
  output_ << "ID=" << header_ptr->getId() << std::endl;
  output_ << "FRACTION_ID=" << header_ptr->getFractionId() << std::endl;
  output_ << "FILE_NAME=" << header_ptr->getFileName() << std::endl;
  output_ << "SCANS=" << header_ptr->getScansString() << std::endl;
  output_ << "RETENTION_TIME=" << std::setprecision(2)
      << header_ptr->getRetentionTime() << std::endl;
  output_ << "LEVEL=" << header_ptr->getMsLevel() << std::endl;

  if (header_ptr->getMsLevel() > 1) {
    if (header_ptr->getActivationPtr() != nullptr) {
      output_ << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
    }
    output_ << "MS_ONE_ID=" << header_ptr->getMsOneId() << std::endl;
    output_ << "MS_ONE_SCAN=" << header_ptr->getMsOneScan() << std::endl;
    output_ << "PRECURSOR_MZ=" << std::setprecision(5) 
        << header_ptr->getPrecMonoMz() << std::endl;
    output_ << "PRECURSOR_CHARGE=" << header_ptr->getPrecCharge() << std::endl;

    /*
    if (use_copied_mono_mass_){//for use in merge sort of msaligns
      output_ << std::fixed << "PRECURSOR_MASS=" << std::setprecision(5) 
        << header_ptr->getCopiedPrecMonoMass() << std::endl;
    }else{
    */
    output_ << std::fixed << "PRECURSOR_MASS=" << std::setprecision(5) 
        << header_ptr->getPrecMonoMass() << std::endl;
    //}
    output_ << "PRECURSOR_INTENSITY=" << std::setprecision(2) 
        <<  header_ptr->getPrecInte() << std::endl;
    /*
    if (header_ptr->getFeatureId() >= 0) {
      output_ << "FEATURE_ID=" << header_ptr->getFeatureId() << std::endl;
      output_ << "FEATURE_INTENSITY=" << std::setprecision(2) 
          << header_ptr->getFeatureInte() << std::endl;
    }
    */
  }

  for (size_t i = 0; i < ms_ptr->size(); i++) {
    DeconvPeakPtr peak_ptr = ms_ptr->getPeakPtr(i);
    output_ << std::setprecision(5) << peak_ptr->getPosition();
    output_ << "\t" << std::setprecision(2) << peak_ptr->getIntensity();
    output_ << "\t" << peak_ptr->getCharge();
    //output_ << "\t" << std::setprecision(2) << peak_ptr->getScore();
    output_ << std::endl;
  }
  output_ << "END IONS" << std::endl;
  output_ << std::endl;
}

}  // namespace toppic
