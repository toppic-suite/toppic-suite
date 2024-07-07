//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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
#include "common/util/str_util.hpp"
#include "common/base/mass_constant.hpp"
#include "common/base/activation_base.hpp"
#include "ms/spec/peak_util.hpp"
#include "ms/spec/msalign_reader.hpp"

namespace toppic {

MsAlignReader::MsAlignReader(const std::string &file_name):
  file_name_(file_name) {
    input_.open(file_name.c_str(), std::ios::in);
    if (!input_.is_open()) {
      LOG_ERROR("msalign file  " << file_name << " does not exist.");
      exit(EXIT_FAILURE);
    }
    group_spec_num_ = 1;
  }

MsAlignReader::MsAlignReader(const std::string &file_name, 
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

MsAlignReader::~MsAlignReader() {
  if (input_.is_open()) {
    input_.close();
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
  int spec_id = -1;

  std::string title;
  std::string scans;
  double retention_time = -1;
  int level = -1;

  int ms_one_id = -1;
  int ms_one_scan = -1;
  double prec_win_begin = -1;
  double prec_win_end = -1;
  std::string activation;

  std::string prec_mass_list = "";
  std::string prec_feat_id_list = "";
  std::string prec_charge_list = "";
  std::string prec_inte_list = "";

  std::vector<std::string> strs;
  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0, 1);
    if (letter >= "A" && letter <= "Z") {
      strs = str_util::split(spectrum_str_vec_[i], "=");
      if (strs[0] == "FILE_NAME") {
        ms_file_name = strs[1];
      } else if (strs[0] == "SPECTRUM_ID") {
        spec_id = std::stoi(strs[1]);
      } else if (strs[0] == "TITLE") {
        title = strs[1];
      } else if (strs[0] == "SCANS") {
        scans = strs[1];
      } else if (strs[0] == "RETENTION_TIME") {
        retention_time = std::stod(strs[1]);
      } else if (strs[0] == "LEVEL") {
        level = std::stoi(strs[1]);
      } else if (strs[0] == "MS_ONE_ID") {
        ms_one_id = std::stoi(strs[1]);
      } else if (strs[0] == "MS_ONE_SCAN") {
        ms_one_scan = std::stod(strs[1]);
      } else if (strs[0] == "PRECURSOR_WINDOW_BEGIN") {
        prec_win_begin = std::stod(strs[1]);
      } else if (strs[0] == "PRECURSOR_WINDOW_END") {
        prec_win_end = std::stod(strs[1]);
      } else if (strs[0] == "ACTIVATION") {
        activation = strs[1];
      } else if (strs[0] == "PRECURSOR_MASS") {
        prec_mass_list = strs[1];
      } else if (strs[0] == "PRECURSOR_CHARGE") {
        prec_charge_list = strs[1];
      } else if (strs[0] == "PRECURSOR_INTENSITY") {
        prec_inte_list = strs[1];
      } else if (strs[0] == "PRECURSOR_FEATURE_ID") {
        prec_feat_id_list = strs[1];
      } 
    }
  }
  MsHeaderPtr header_ptr = std::make_shared<MsHeader>();
  header_ptr->setFileName(ms_file_name);
  //set spec_id
  if (spec_id < 0) {
    LOG_ERROR("Spectrum id is missing!");
    exit(EXIT_FAILURE);
  }
  header_ptr->setSpecId(spec_id);
  // set title
  if (title == "") {
    title = "sp_" + str_util::toString(spec_id);
  }
  header_ptr->setTitle(title);
  // set scans
  if (scans != "") {
    header_ptr->setScans(scans);
  } else {
    header_ptr->setScans("");
  }
  header_ptr->setRetentionTime(retention_time);
  // set ms level
  if (level <= 0) {
    LOG_ERROR("MS level information is missing in MSALIGN file!");
    exit(EXIT_FAILURE);
  }
  header_ptr->setMsLevel(level);

  if (level > 1) {
    header_ptr->setMsOneId(ms_one_id);
    header_ptr->setMsOneScan(ms_one_scan);
    //set prec window 
    if (prec_win_begin < 0 || prec_win_end < 0) {
      LOG_ERROR("Precursor window information is missing in MSALIGN file!");
      exit(EXIT_FAILURE);
    }
    header_ptr->setPrecWinBegin(prec_win_begin);
    header_ptr->setPrecWinEnd(prec_win_end);
    // set activation type
    if (activation_ptr_ != nullptr) {
      // use the default activation if the information is missing
      header_ptr->setActivationPtr(activation_ptr_);
    } else if (activation != "") {
      ActivationPtr activation_ptr = ActivationBase::getActivationPtrByName(activation);
      header_ptr->setActivationPtr(activation_ptr);
    }
    // set precursor information
    PrecursorPtrVec prec_ptr_vec;
    if (prec_mass_list != "") {
      std::vector<std::string> mass_strs = str_util::split(prec_mass_list, ":");
      std::vector<std::string> feat_id_strs = str_util::split(prec_feat_id_list, ":");
      std::vector<std::string> charge_strs = str_util::split(prec_charge_list, ":");
      std::vector<std::string> inte_strs = str_util::split(prec_inte_list, ":");
      for (size_t id = 0; id < mass_strs.size(); id++) {    
        double prec_mass = std::stod(mass_strs[id]);
        int prec_feat_id = std::stoi(feat_id_strs[id]);
        int prec_charge = std::stoi(charge_strs[id]);
        double prec_inte = std::stod(inte_strs[id]);
        double prec_mono_mz = peak_util::compMz(prec_mass, prec_charge); 
        PrecursorPtr prec_ptr = std::make_shared<Precursor>(id, prec_feat_id, 
                                                            prec_mono_mz, prec_charge, prec_inte);
        prec_ptr_vec.push_back(prec_ptr);
      }
    }
    header_ptr->setPrecPtrVec(prec_ptr_vec);
  }

  std::vector<DeconvPeakPtr> peak_ptr_list;
  int peak_id = 0;
  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0, 1);
    if (letter >= "0" && letter <= "9") {
      strs = str_util::split(spectrum_str_vec_[i], "\t ");
      double mass = std::stod(strs[0]);
      double inte = std::stod(strs[1]);
      double charge = std::stoi(strs[2]);
      double score = 1.0;
      if (strs.size() > 3) {
        score = std::stod(strs[3]);
      }
      DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(spec_id, peak_id, mass, inte,
                                                            charge, score);
      peak_ptr_list.push_back(peak_ptr);
      peak_id++;
    }
  }

  deconv_ms_ptr_ = std::make_shared<Ms<DeconvPeakPtr> >(header_ptr, peak_ptr_list);
  current_++;
}

DeconvMsPtr MsAlignReader::getNextMsPtr() {
  readNext();
  return deconv_ms_ptr_;
}

DeconvMsPtrVec MsAlignReader::getNextMsPtrVec() {
  DeconvMsPtrVec deconv_ms_ptr_vec;
  for (int i = 0; i < group_spec_num_; i++) {
    readNext();
    if (deconv_ms_ptr_ == nullptr) {
      deconv_ms_ptr_vec.clear();
      return deconv_ms_ptr_vec;
    }
    deconv_ms_ptr_vec.push_back(deconv_ms_ptr_);
  }
  // make sure that all MS/MS spectra have the same precursors 
  // we may need to implement deep copy here to make sure the 
  // scans do not share precursor instances
  PrecursorPtrVec prec_ptr_vec =
    deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecPtrVec();
  for (int i = 1; i < group_spec_num_; i++) {
    deconv_ms_ptr_vec[i]->getMsHeaderPtr()->setPrecPtrVec(prec_ptr_vec);
  }
  return deconv_ms_ptr_vec;
}

}  // namespace toppic
