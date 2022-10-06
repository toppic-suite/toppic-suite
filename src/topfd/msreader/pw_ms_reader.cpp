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

#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/data/common/cv.hpp"
#include "pwiz/data/common/ParamTypes.hpp"

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/base/activation_base.hpp"
#include "topfd/msreader/pw_ms_reader.hpp"

namespace toppic {
PwMsReader::PwMsReader(const std::string & file_name,
                       double isolation_window):
    file_name_(file_name),
    isolation_window_(isolation_window),
    input_sp_id_(0),
    output_sp_id_(0) {
      msd_ptr_ = std::make_shared<pwiz::msdata::MSDataFile>(file_name, &readers_);
      spec_list_ptr_ =  msd_ptr_->run.spectrumListPtr;
      input_sp_num_ = spec_list_ptr_->size();
    }

PwMsReader::PwMsReader(const std::string & file_name, std::string activation,
                       double isolation_window):
    file_name_(file_name),
    activation_(activation),
    isolation_window_(isolation_window),
    input_sp_id_(0),
    output_sp_id_(0) {
      msd_ptr_ = std::make_shared<pwiz::msdata::MSDataFile>(file_name, &readers_);
      spec_list_ptr_ =  msd_ptr_->run.spectrumListPtr;
      input_sp_num_ = spec_list_ptr_->size();
    }


bool PwMsReader::readOneMs(int sp_id, PeakPtrVec &peak_list, MsHeaderPtr &header_ptr) {
  pwiz::msdata::SpectrumPtr cur_spec_ptr = nullptr;
  bool get_binary_data = true;
  cur_spec_ptr = spec_list_ptr_->spectrum(sp_id, get_binary_data);
  if (cur_spec_ptr == nullptr) {return false;}

  bool is_centroided = cur_spec_ptr->hasCVParam(pwiz::cv::MS_centroid_spectrum);
  if (!is_centroided) {
    std::cout << "Error: The data file contains profile data!" << std::endl;
    std::cout << "TopFD can process only centroided, not profile, MS data!" << std::endl;  
    exit(EXIT_FAILURE);
  }

  // get m/z and intensity values from the spectra
  std::vector<pwiz::msdata::MZIntensityPair> pairs;
  cur_spec_ptr->getMZIntensityPairs(pairs);
  LOG_DEBUG("mz pair size " << pairs.size());
  // make sure the peak list is already sorted
  std::sort(pairs.begin(), pairs.end(),
            [](const pwiz::msdata::MZIntensityPair &a, const pwiz::msdata::MZIntensityPair& b) {
            return a.mz < b.mz;
            });
  peak_list.empty();
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i].intensity > 0.0) {
      PeakPtr peak_ptr = std::make_shared<Peak>(pairs[i].mz, pairs[i].intensity);
      peak_list.push_back(peak_ptr);
    }
  }

  // scan info
  pwiz::msdata::SpectrumInfo spec_info(*cur_spec_ptr);
  int ms_level = spec_info.msLevel;
  LOG_DEBUG("ms_level " << ms_level);
  // For agilent data, scan numbers are missing.
  if (spec_info.scanNumber == 0) {
    spec_info.scanNumber = spec_info.index + 1;
  }
  LOG_DEBUG(" scan number " << spec_info.scanNumber);

  if (ms_level == 2) {
    header_ptr = std::make_shared<MsHeader>();
    header_ptr->setId(ms2_cnt);
    ms2_cnt++;
    header_ptr->setScan(spec_info.scanNumber);
    header_ptr->setMsLevel(ms_level);
    header_ptr->setFileName(file_name_);
    header_ptr->setTitle("Scan_" + str_util::toString(spec_info.scanNumber));
    header_ptr->setRetentionTime(spec_info.retentionTime);

    // precursor information
    double prec_mz = 0;
    int prec_charge = 1;
    double prec_inte = 0.0;
    double prec_scan_num = -1;
    if (spec_info.precursors.size() > 0) {
      prec_mz = spec_info.precursors[0].mz;
      prec_charge = static_cast<int>(spec_info.precursors[0].charge);
      prec_inte = spec_info.precursors[0].intensity;
      if (prec_mz < 0) {prec_mz = 0;}
      if (prec_charge  < 0) {prec_charge = 1;}
      if (prec_inte < 0) {prec_inte = 0.0;}
      //get precursor scan ID from mzML
      std::string line = cur_spec_ptr->precursors[0].spectrumID;
      std::string delimiter = "=";
      if (line.substr(line.find_last_of(delimiter) + 1) == "") {
        prec_scan_num = prev_ms1_scan_id;
      }
      else {
        prec_scan_num = std::stoi(line.substr(line.find_last_of(delimiter) + 1)); 
      }
    }
    // precursor mz in mzML data
    header_ptr->setPrecSpMz(prec_mz);
    header_ptr->setPrecCharge(prec_charge);
    header_ptr->setPrecInte(prec_inte);
    header_ptr->setMsOneScan(prec_scan_num);
    LOG_DEBUG("Precursor m/z " << prec_mz);

    // isolation window
    double prec_target_mz = prec_mz;
    double isolation_lower_offset = isolation_window_ / 2;
    double isolation_upper_offset = isolation_window_/ 2;
    if (cur_spec_ptr->precursors.size() > 0) {
      std::vector<pwiz::data::CVParam> cv_list = cur_spec_ptr->precursors[0].isolationWindow.cvParams;
      for (size_t i = 0; i < cv_list.size(); i++) {
        LOG_DEBUG("cv list " << i << " " << cv_list[i].cvid << " " << cv_list[i].value);
        if (cv_list[i].cvid == pwiz::cv::MS_isolation_window_target_m_z) {
          prec_target_mz = std::stod(cv_list[i].value);
        }
        if (cv_list[i].cvid == pwiz::cv::MS_isolation_window_lower_offset) {
          isolation_lower_offset = std::stod(cv_list[i].value);
        }
        if (cv_list[i].cvid == pwiz::cv::MS_isolation_window_upper_offset) {
          isolation_upper_offset = std::stod(cv_list[i].value);
        }
      }
    }
    LOG_DEBUG("Precursor target mz " << prec_target_mz << " lower offset " << isolation_lower_offset 
              << " upper offset " << isolation_upper_offset); 
    header_ptr->setPrecTargetMz(prec_target_mz);
    header_ptr->setIsolationLowerOffset(isolation_lower_offset);
    header_ptr->setIsolationUpperOffset(isolation_upper_offset);

    // activation type
    std::string ac_name = activation_;
    if (ac_name == "" || ac_name == "FILE"){
      ac_name = "";
      if (cur_spec_ptr->precursors.size() > 0) {
        std::vector<pwiz::data::CVParam> cv_list = cur_spec_ptr->precursors[0].activation.cvParams;
        for (size_t i = 0; i < cv_list.size(); i++) {
          LOG_DEBUG("cv list " << i << " " << cv_list[i].cvid << " " << cv_list[i].value);
          if (cv_list[i].cvid == pwiz::cv::MS_CID) {
            ac_name = "CID";
            break;
          } else if (cv_list[i].cvid == pwiz::cv::MS_HCD) {
            ac_name = "HCD";
            break;
          } else if (cv_list[i].cvid == pwiz::cv::MS_ETD) {
            ac_name = "ETD";
            break;
          } else if (cv_list[i].cvid == pwiz::cv::MS_MPD) {
            ac_name = "MPD";
            break;
          }
        }
      }
    }
    if (ac_name == "") {
      LOG_WARN("No activation information is available in reading the spectrum with scan " << spec_info.scanNumber);
      std::cout << "\nERROR: Unable to read the activation method from the file.";
      std::cout << "\nPlease select an activation method out of CID|HCD|ETD|UVPD and provide it as a parameter.";
      std::cout << "\nExample: -a CID" << std::endl;
      exit(EXIT_FAILURE);
    }
    LOG_DEBUG("ac name " << ac_name);
    ActivationPtr activation_ptr = ActivationBase::getActivationPtrByName(ac_name);
    header_ptr->setActivationPtr(activation_ptr);
  } 
  else {
    header_ptr = std::make_shared<MsHeader>();
    header_ptr->setId(ms1_cnt);
    ms1_cnt++;
    header_ptr->setScan(spec_info.scanNumber);
    header_ptr->setMsLevel(ms_level);
    header_ptr->setPrecCharge(0);
    header_ptr->setFileName(file_name_);
    header_ptr->setTitle("Scan_" + str_util::toString(spec_info.scanNumber));
    header_ptr->setRetentionTime(spec_info.retentionTime);
    prev_ms1_scan_id = spec_info.scanNumber;
  }
  //add voltage information if it exists
  std::vector<pwiz::data::CVParam> cv_list = (*cur_spec_ptr).cvParams;
  for (size_t i = 0; i < cv_list.size(); i++) {
    if (cv_list[i].cvid == pwiz::cv::MS_FAIMS_compensation_voltage) {
      header_ptr->setVoltage(std::stod(cv_list[i].value));
      break;
    }
  }
  return true;
}

int PwMsReader::readNext() {
  peak_list_.clear();
  header_ptr_ = nullptr;

  int found = false;
  while (!found) {
    if (input_sp_id_ >= input_sp_num_) {
      LOG_DEBUG("Only " << input_sp_num_  << " spectra in the input data!");
      return -1;
    }
    found = readOneMs(input_sp_id_, peak_list_, header_ptr_); 
    input_sp_id_++;
  }
  return 1;
}

}  // namespace toppic
