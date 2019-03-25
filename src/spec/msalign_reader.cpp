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

#include <cmath>
#include <sstream>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/base/mass_constant.hpp"
#include "common/base/activation_base.hpp"
#include "spec/msalign_reader.hpp"

namespace toppic {


MsAlignReader::MsAlignReader(const std::string &file_name):
    file_name_(file_name) {
      input_.open(file_name.c_str(), std::ios::in);
      if (!input_.is_open()) {
        LOG_ERROR("msalign file  " << file_name << " does not exist.");
        exit(EXIT_FAILURE);
      }
      group_spec_num_ = 1;
      activation_ptr_ = nullptr;
      peak_num_limit_ = std::numeric_limits<int>::max();
    }
      
MsAlignReader::MsAlignReader(const std::string &file_name, int group_spec_num,
                             ActivationPtr act_ptr, const std::set<std::string> skip_list,
                             int peak_num_limit):
    file_name_(file_name),
    group_spec_num_(group_spec_num),
    activation_ptr_(act_ptr),
    skip_list_(skip_list),
    peak_num_limit_(peak_num_limit) {
      input_.open(file_name.c_str(), std::ios::in);
      if (!input_.is_open()) {
        LOG_ERROR("msalign file  " << file_name << " does not exist.");
        exit(EXIT_FAILURE);
      }
    }

std::vector<std::string> MsAlignReader::readOneStrSpectrum() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    str_util::trim(line);
    if (line == "BEGIN IONS") {
      line_list.push_back(line);
    } else if (line == "END IONS") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    } else if (line == "" || line[0] == '#') {
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
  deconv_ms_ptr_ = nullptr;
  spectrum_str_vec_ = readOneStrSpectrum();
  if (spectrum_str_vec_.size() == 0) {
    input_.close();
    return;
  }
  std::string ms_file_name = "";
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
    std::string letter = spectrum_str_vec_[i].substr(0, 1);
    if (letter >= "A" && letter <= "Z") {
      strs = str_util::split(spectrum_str_vec_[i], "=");
      if (strs[0] == "ID") {
        id = std::stoi(strs[1]);
      }

      if (strs[0] == "FILE_NAME") {
        ms_file_name = strs[1];
      } else if (strs[0] == "PRECURSOR_ID") {
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

  MsHeaderPtr header_ptr = std::make_shared<MsHeader>();
  header_ptr->setFileName(ms_file_name);
  header_ptr->setId(id);
  header_ptr->setPrecId(prec_id);

  if (scans != "") {
    header_ptr->setScans(scans);
  } else {
    header_ptr->setScans("");
  }
  header_ptr->setRetentionTime(retention_time);
  // LOG_DEBUG("retention time " << retention_time);

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
    ActivationPtr activation_ptr = ActivationBase::getActivationPtrByName(activation);

    header_ptr->setActivationPtr(activation_ptr);
  }
  header_ptr->setMsLevel(level);

  header_ptr->setMsOneId(ms_one_id);

  header_ptr->setMsOneScan(ms_one_scan);

  header_ptr->setPrecMonoMz(prec_mass /prec_charge + mass_constant::getProtonMass());

  header_ptr->setPrecCharge(prec_charge);

  header_ptr->setPrecInte(prec_inte);

  header_ptr->setFeatureId(feature_id);

  header_ptr->setFeatureInte(feature_inte);

  std::vector<DeconvPeakPtr> peak_ptr_list;
  int idx = 0;
  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0, 1);
    if (letter >= "0" && letter <= "9") {
      strs = str_util::split(spectrum_str_vec_[i], "\t ");
      double mass = std::stod(strs[0]);
      double inte = std::stod(strs[1]);
      int charge = std::stoi(strs[2]);
      DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(id, idx, mass, inte, charge);
      peak_ptr_list.push_back(peak_ptr);
      idx++;
      if (static_cast<int>(peak_ptr_list.size()) >= peak_num_limit_) break;
    }
  }

  deconv_ms_ptr_ = std::make_shared<Ms<DeconvPeakPtr> >(header_ptr, peak_ptr_list);

  current_++;

}

DeconvMsPtr MsAlignReader::getNextMs() {
  readNext();
  while (deconv_ms_ptr_ != nullptr
         && skip_list_.find(deconv_ms_ptr_->getMsHeaderPtr()->getScansString()) != skip_list_.end()) {
    readNext();
  }
  return deconv_ms_ptr_;
}

std::vector<SpectrumSetPtr> MsAlignReader::getNextSpectrumSet(SpParaPtr sp_para_ptr) {
  std::vector<SpectrumSetPtr> spec_set_vec;
  DeconvMsPtrVec deconv_ms_ptr_vec;
  for (int i = 0; i < group_spec_num_; i++) {
    readNext();
    while (deconv_ms_ptr_ != nullptr
           && skip_list_.find(deconv_ms_ptr_->getMsHeaderPtr()->getScansString()) != skip_list_.end()) {
      readNext();
    }
    if (deconv_ms_ptr_ == nullptr) {
      spec_set_vec.push_back(nullptr);
      return spec_set_vec;
    }
    deconv_ms_ptr_vec.push_back(deconv_ms_ptr_);
  }
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  // LOG_DEBUG("prec mass " << prec_mono_mass);
  int count = 1;
  for (int i = 1; i < group_spec_num_; i++) {
    double new_mass = deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getPrecMonoMass();
    if (std::abs(prec_mono_mass - new_mass) < 0.5) {
      prec_mono_mass = (prec_mono_mass * count + new_mass)/ (count+1);
      count++;
    }
  }
  // LOG_DEBUG("prec mass result " << prec_mono_mass);
  std::vector<double> prec_errors;
  prec_errors.push_back(0);
  for (int i = 1; i <= sp_para_ptr->prec_error_; i++) {
    prec_errors.push_back(- i * mass_constant::getIsotopeMass());
    prec_errors.push_back(i * mass_constant::getIsotopeMass());
  }
  for (size_t i = 0; i< prec_errors.size(); i++) {
    spec_set_vec.push_back(std::make_shared<SpectrumSet>(deconv_ms_ptr_vec,
                                                         sp_para_ptr,
                                                         prec_mono_mass + prec_errors[i]));
  }
  return spec_set_vec;
}

void MsAlignReader::close() {
  input_.close();
}

}  // namespace toppic
