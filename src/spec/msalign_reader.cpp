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


#include <cmath>

#include "base/logger.hpp"
#include "base/activation_base.hpp"
#include "base/string_util.hpp"
#include "spec/msalign_reader.hpp"

namespace prot {

MsAlignReader::MsAlignReader(const std::string &file_name, 
                             int group_spec_num, ActivationPtr act_ptr): 
    file_name_(file_name),
    group_spec_num_(group_spec_num),
    activation_ptr_(act_ptr) {
      input_.open(file_name.c_str(), std::ios::in);
    }

std::vector<std::string> MsAlignReader::readOneSpectrum() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    line = StringUtil::trim(line);
    if (line ==  "BEGIN IONS") {
      line_list.push_back(line);
    } else if (line == "END IONS") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    } else if (line == "") {
      continue;
    } else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

void MsAlignReader::readNext() {
  deconv_ms_ptr_ = DeconvMsPtr(nullptr);
  spectrum_str_vec_ = readOneSpectrum();
  if (spectrum_str_vec_.size() == 0) {
    input_.close();
    return;
  }
  int id = -1;
  int prec_id = 0;
  std::string scans;
  double retention_time = -1;
  std::string activation;
  std::string title;
  int level = 2;
  int ms_one_id = -1;
  int ms_one_scan = -1;
  double prec_mass = -1;
  int prec_charge = -1;
  double prec_inte = -1;
  int feature_id = -1;
  double feature_inte = -1;
  std::vector<std::string> strs;
  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0,1);
    if (letter >= "A" && letter <= "Z") {
      strs = StringUtil::split(spectrum_str_vec_[i], '=');
      if (strs[0] == "ID") {
        id = std::stoi(strs[1]);
      }

      if (strs[0] == "PRECURSOR_ID") {
        prec_id = std::stoi(strs[1]);
      } else if (strs[0] == "SCANS") {
        scans = strs[1];
      } else if (strs[0] == "RETENTION_TIME") {
        retention_time = std::stod(strs[1]);        
      } else if (strs[0] == "ACTIVATION") {
        activation = strs[1];
      } else if (strs[0] == "TITLE") {
        title = strs[1];
      } else if (strs[0] == "LEVEL") {
        level = std::stoi(strs[1]);
      } else if (strs[0] == "MS_ONE_ID") {
        ms_one_id = std::stod(strs[1]);
      } else if (strs[0] == "MS_ONE_SCAN") {
        ms_one_scan = std::stod(strs[1]);
      } else if (strs[0] == "PRECURSOR_MASS") {
        prec_mass = std::stod(strs[1]);
      } else if (strs[0] == "PRECURSOR_CHARGE") {
        prec_charge = std::stoi(strs[1]);
      } else if (strs[0] == "PRECURSOR_INTENSITY") {
        prec_inte = std::stod(strs[1]);
      } else if (strs[0] == "FEATURE_ID") {
        feature_id = std::stoi(strs[1]);
      } else if (strs[0] == "FEATURE_INTENSITY") {
        feature_inte = std::stod(strs[1]);
      }
    }
  }
  if (id < 0 || prec_charge < 0 || prec_mass < 0) {
    LOG_WARN("Input file format error: sp id " << id << " prec_chrg "
              << prec_charge << " prec mass " << prec_mass);
  }

  MsHeaderPtr header_ptr(new MsHeader());
  header_ptr->setFileName(file_name_);
  header_ptr->setId(id);
  header_ptr->setPrecId(prec_id);
  if (scans != "") {
    header_ptr->setScans(scans);
  } else {
    header_ptr->setScans("");
  }
  header_ptr->setRetentionTime(retention_time);
  //LOG_DEBUG("retention time " << retention_time);

  if (title != "") {
    std::stringstream ss;
    ss << "sp_" << id;
    header_ptr->setTitle(ss.str());
  } else {
    header_ptr->setTitle(title);
  }

  if (activation_ptr_ != nullptr) {
    header_ptr->setActivationPtr(activation_ptr_);
  } else if (activation != "") {
    ActivationPtr activation_ptr = 
        ActivationBase::getActivationPtrByName(activation);
    header_ptr->setActivationPtr(activation_ptr);
  }
  header_ptr->setMsLevel(level);

  header_ptr->setMsOneId(ms_one_id);

  header_ptr->setMsOneScan(ms_one_scan);

  header_ptr->setPrecMonoMz(prec_mass /prec_charge
                            + MassConstant::getProtonMass());
  header_ptr->setPrecCharge(prec_charge);
  
  header_ptr->setPrecInte(prec_inte);

  header_ptr->setFeatureId(feature_id);

  header_ptr->setFeatureInte(feature_inte);

  std::vector<DeconvPeakPtr> peak_ptr_list;
  int idx = 0;
  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0,1);
    if (letter >= "0" && letter <= "9") {
      strs = StringUtil::split(spectrum_str_vec_[i], '\t');
      double mass = std::stod(strs[0]);
      double inte = std::stod(strs[1]);
      int charge = std::stoi(strs[2]);
      DeconvPeakPtr peak_ptr(new DeconvPeak(idx, mass, inte, charge));
      peak_ptr_list.push_back(peak_ptr);
      idx++;
    }
  }
  deconv_ms_ptr_ 
      = DeconvMsPtr(new Ms<DeconvPeakPtr>(header_ptr, peak_ptr_list));

  current_++;
}

DeconvMsPtr MsAlignReader::getNextMs() {
  readNext();
  return deconv_ms_ptr_;
}

SpectrumSetPtr MsAlignReader::getNextSpectrumSet(SpParaPtr sp_para_ptr) {
  DeconvMsPtrVec deconv_ms_ptr_vec; 
  for (int i = 0; i < group_spec_num_; i++) {
    readNext();
    if (deconv_ms_ptr_ == nullptr) {
      return SpectrumSetPtr(nullptr);
    }
    deconv_ms_ptr_vec.push_back(deconv_ms_ptr_);
  }
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  //LOG_DEBUG("prec mass " << prec_mono_mass);
  int count = 1;
  for (int i = 1; i < group_spec_num_; i++) {
    double new_mass = deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getPrecMonoMass();
    if (std::abs(prec_mono_mass - new_mass) < 0.5) {
      prec_mono_mass = (prec_mono_mass * count + new_mass)/ (count+1);
      count++;
    }
  }
  //LOG_DEBUG("prec mass result " << prec_mono_mass);
  return SpectrumSetPtr(new SpectrumSet(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass));
}

void MsAlignReader::close() {
  input_.close();
}

}
