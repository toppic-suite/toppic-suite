//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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
    LOG_ERROR("Cannot open the msalign file  " << file_name << "!");
    exit(EXIT_FAILURE);
  }
}

MsAlignWriter::~MsAlignWriter() {
  if (output_.is_open()) {
    output_.close();
  }
}

void MsAlignWriter::writePara(const std::string &para_str) {
  output_ << para_str << "\n";
}

void MsAlignWriter::writeMs(DeconvMsPtr ms_ptr) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  output_ << std::fixed;
  output_ << "BEGIN IONS" << std::endl;
  output_ << "FILE_NAME=" << header_ptr->getFileName() << std::endl;
  output_ << "SPECTRUM_ID=" << header_ptr->getSpecId() << std::endl;
  output_ << "TITLE=" << header_ptr->getTitle() << std::endl;
  output_ << "SCANS=" << header_ptr->getScansString() << std::endl;
  output_ << "RETENTION_TIME=" << std::fixed << std::setprecision(2)
      << header_ptr->getRetentionTime() << std::endl;
  output_ << "LEVEL=" << header_ptr->getMsLevel() << std::endl;

  if (header_ptr->getMsLevel() > 1) {
    output_ << "MS_ONE_ID=" << header_ptr->getMsOneId() << std::endl;
    output_ << "MS_ONE_SCAN=" << header_ptr->getMsOneScan() << std::endl;
    output_ << "PRECURSOR_WINDOW_BEGIN=" << header_ptr->getPrecWinBegin() << std::endl;
    output_ << "PRECURSOR_WINDOW_END=" << header_ptr->getPrecWinEnd() << std::endl;
    if (header_ptr->getActivationPtr() != nullptr) {
      output_ << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
    }
    PrecursorPtrVec prec_ptrs = header_ptr->getPrecPtrVec();
    output_ << "PRECURSOR_MZ=" << std::fixed << std::setprecision(5);
    for (size_t i = 0; i < prec_ptrs.size(); i++) {
      output_ << prec_ptrs[i]->getMonoMz();
      // use : for separating multiple precursors
      if (i < prec_ptrs.size() - 1) {output_ << ":";}
    }
    output_ << std::endl;
    output_ << "PRECURSOR_CHARGE=";
    for (size_t i = 0; i < prec_ptrs.size(); i++) {
      output_ << prec_ptrs[i]->getCharge();
      // use : for separating multiple precursors
      if (i < prec_ptrs.size() - 1) {output_ << ":";}
    }
    output_ << std::endl;
    // The precision for mass is 5
    output_ << "PRECURSOR_MASS=" << std::fixed << std::setprecision(5);
    for (size_t i = 0; i < prec_ptrs.size(); i++) {
      output_ << prec_ptrs[i]->getMonoMass();
      // use : for separating multiple precursors
      if (i < prec_ptrs.size() - 1) {output_ << ":";}
    }
    output_ << std::endl;
    output_ << "PRECURSOR_INTENSITY=" << std::fixed << std::setprecision(2);
    for (size_t i = 0; i < prec_ptrs.size(); i++) {
      output_ << prec_ptrs[i]->getInte();
      // use : for separating multiple precursors
      if (i < prec_ptrs.size() - 1) {output_ << ":";}
    }
    output_ << std::endl;
    output_ << "PRECURSOR_FEATURE_ID=";
    for (size_t i = 0; i < prec_ptrs.size(); i++) {
      output_ << prec_ptrs[i]->getFeatureId();
      // use : for separating multiple precursors
      if (i < prec_ptrs.size() - 1) {output_ << ":";}
    }
    output_ << std::endl;
  }
  for (size_t i = 0; i < ms_ptr->size(); i++) {
    DeconvPeakPtr peak_ptr = ms_ptr->getPeakPtr(i);
    output_ << std::fixed << std::setprecision(5) << peak_ptr->getPosition();
    output_ << "\t" << std::fixed << std::setprecision(2) << peak_ptr->getIntensity();
    output_ << "\t" << peak_ptr->getCharge();
    output_ << "\t" << std::fixed << std::setprecision(2) << peak_ptr->getScore();
    output_ << std::endl;
  }
  output_ << "END IONS" << std::endl;
  output_ << std::endl;
}

}  // namespace toppic
