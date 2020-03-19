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

#include <cmath>
#include <sstream>
#include <iomanip>
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/base/mass_constant.hpp"
#include "common/base/activation_base.hpp"
#include "ms/spec/simple_msalign_reader.hpp"

namespace toppic {

SimpleMsAlignReader::SimpleMsAlignReader(const std::string &file_name):
    file_name_(file_name) {
      input_.open(file_name.c_str(), std::ios::in);
      if (!input_.is_open()) {
        LOG_ERROR("msalign file  " << file_name << " does not exist.");
        exit(EXIT_FAILURE);
      }
      group_spec_num_ = 1;
    }

SimpleMsAlignReader::SimpleMsAlignReader(const std::string &file_name, 
                                         int group_spec_num,
                                         ActivationPtr activation_ptr):
    file_name_(file_name),
    group_spec_num_(group_spec_num),
    activation_ptr_(activation_ptr) {
      input_.open(file_name.c_str(), std::ios::in);
      if (!input_.is_open()) {
        LOG_ERROR("msalign file  " << file_name << " does not exist.");
        exit(EXIT_FAILURE);
      }
    }

SimpleMsAlignReader::~SimpleMsAlignReader() {
  if (input_.is_open()) {
    input_.close();
  }
}
      
std::vector<std::string> SimpleMsAlignReader::readOneStrSpectrum() {
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

void SimpleMsAlignReader::readNext() {
  deconv_ms_ptr_ = nullptr;
  spectrum_str_vec_ = readOneStrSpectrum();
  if (spectrum_str_vec_.size() == 0) {
    input_.close();
    return;
  }
  std::string ms_file_name = "";
  int id = -1;
  int fraction_id = -1;
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
  //int feature_id = -1;
  //double feature_inte = -1;
  std::vector<std::string> strs;

  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0, 1);
    if (letter >= "A" && letter <= "Z") {
      strs = str_util::split(spectrum_str_vec_[i], "=");
      if (strs[0] == "ID") {
        id = std::stoi(strs[1]);
      } else if (strs[0] == "FRACTION_ID") {
        fraction_id = std::stoi(strs[1]);
      } else if (strs[0] == "FILE_NAME") {
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
      } 
      /*
      else if (strs[0] == "FEATURE_ID") {
        feature_id = std::stoi(strs[1]);
      } else if (strs[0] == "FEATURE_INTENSITY") {
        feature_inte = std::stod(strs[1]);
      }
      */
    }
  }
  if (id < 0 || prec_charge < 0 || prec_mass < 0) {
    if (level == 2) {
      LOG_WARN("Input file format error: sp id " << id << " prec_chrg "
               << prec_charge << " prec mass " << prec_mass);
    }
  }

  MsHeaderPtr header_ptr = std::make_shared<MsHeader>();
  header_ptr->setFractionId(fraction_id);
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

  //header_ptr->setFeatureId(feature_id);

  //header_ptr->setFeatureInte(feature_inte);

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
    }
  }

  deconv_ms_ptr_ = std::make_shared<Ms<DeconvPeakPtr> >(header_ptr, peak_ptr_list);
  current_++;
}

DeconvMsPtr SimpleMsAlignReader::getNextMsPtr() {
  readNext();
  return deconv_ms_ptr_;
}

DeconvMsPtrVec SimpleMsAlignReader::getNextMsPtrVec() {
  DeconvMsPtrVec deconv_ms_ptr_vec;
  for (int i = 0; i < group_spec_num_; i++) {
    readNext();
    if (deconv_ms_ptr_ == nullptr) {
      deconv_ms_ptr_vec.clear();
      return deconv_ms_ptr_vec;
    }
    deconv_ms_ptr_vec.push_back(deconv_ms_ptr_);
  }
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  int charge = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge();
  // Get average precursor mass
  int count = 1;
  for (int i = 1; i < group_spec_num_; i++) {
    double new_mass = deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getPrecMonoMass();
    if (std::abs(prec_mono_mass - new_mass) < 0.5) {
      prec_mono_mass = (prec_mono_mass * count + new_mass)/ (count+1);
      count++;
    }
  }

  double prec_mono_mz = Peak::compMz(prec_mono_mass, charge);
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    deconv_ms_ptr_vec[i]->getMsHeaderPtr()->setPrecMonoMz(prec_mono_mz);
  }
  return deconv_ms_ptr_vec;
}

void SimpleMsAlignReader::readMsOneSpectra(const std::string &file_name, 
                                           DeconvMsPtrVec &ms_ptr_vec) {
  SimpleMsAlignReader sp_reader(file_name); 

  DeconvMsPtr ms_ptr;
  while ((ms_ptr = sp_reader.getNextMsPtr())!= nullptr) {
    ms_ptr->getMsHeaderPtr()->setMsLevel(1);
    ms_ptr_vec.push_back(ms_ptr);
  }
}

}  // namespace toppic
