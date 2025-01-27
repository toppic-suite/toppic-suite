//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include <set>

#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/data/common/cv.hpp"
#include "pwiz/data/common/ParamTypes.hpp"

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/base/activation_base.hpp"
#include "ms/mzml/pw_ms_reader.hpp"

namespace toppic {

PwMsReader::PwMsReader(const std::string & file_name) {
  init(file_name);
}

PwMsReader::PwMsReader(const std::string & file_name,
                       double isolation_window) {
  init(file_name);
  isolation_window_ = isolation_window;
}

PwMsReader::PwMsReader(const std::string & file_name, 
                       double isolation_window, 
                       std::string activation) {
  init(file_name);
  isolation_window_ = isolation_window;
  activation_ = activation;
}

void PwMsReader::init(const std::string & file_name) {
  file_name_ = file_name;
  input_sp_id_ = 0;
  msd_ptr_ = std::make_shared<pwiz::msdata::MSDataFile>(file_name, &readers_);
  spec_list_ptr_ =  msd_ptr_->run.spectrumListPtr;
  input_sp_num_ = spec_list_ptr_->size();
  is_waters_instrument_ = checkWatersInstrument();
}

bool PwMsReader::checkWatersInstrument() {
  for (size_t i = 0; i < msd_ptr_->instrumentConfigurationPtrs.size(); i++) {
    pwiz::msdata::InstrumentConfigurationPtr conf_ptr = msd_ptr_->instrumentConfigurationPtrs[i];
    pwiz::msdata::CVParam param = conf_ptr->cvParam(pwiz::msdata::CVID::MS_Waters_instrument_model);
    if (! param.empty()) {
      LOG_DEBUG("Found waters model"); 
      return true; 
    }
  }
  return false;
}

bool PwMsReader::checkCentroidData() {
  std::vector<pwiz::msdata::DataProcessingPtr> proc_ptr_vec = msd_ptr_->dataProcessingPtrs;
  for (size_t i = 0; i < proc_ptr_vec.size(); i++) {
    for (size_t j = 0; j < proc_ptr_vec[i]->processingMethods.size(); j++) {
      pwiz::msdata::ProcessingMethod proc_method = proc_ptr_vec[i]->processingMethods[j];
      pwiz::cv::CVID cvid_peak_picking;
      cvid_peak_picking = pwiz::cv::MS_peak_picking;
      if (proc_method.hasCVParam(cvid_peak_picking)) {
        return true;
      }
    }
  }
  return false;
}

void PwMsReader::resetIndexes() {
  input_sp_id_ = 0;
  ms1_cnt_ = 0;
  ms2_cnt_ = 0;
  prev_ms1_scan_id_ = -1;
}

int PwMsReader::parseNum(std::string &id, int default_scan) {
  std::string delimiter = "=";
  // 
  if (id.substr(id.find_last_of(delimiter) + 1) == "") {
    LOG_INFO("Scan number information is missing!");
    return default_scan;
  }
  else {
    return std::stoi(id.substr(id.find_last_of(delimiter) + 1)); 
  }
}

void PwMsReader::parseScanNum(MsHeaderPtr header_ptr, 
                              pwiz::msdata::SpectrumInfo &spec_info) {
  // For Agilent data, scan numbers are missing.
  if (spec_info.scanNumber == 0) {
    // default is index 
    spec_info.scanNumber = parseNum(spec_info.id, spec_info.index);
  }

  // For Waters data, use spec_info.index as the scan number
  if (is_waters_instrument_) {
    spec_info.scanNumber = spec_info.index;
  }
  LOG_DEBUG(" scan number " << spec_info.scanNumber);
  header_ptr->setSingleScan(spec_info.scanNumber);
}

PeakPtrVec PwMsReader::parsePeaks(pwiz::msdata::SpectrumPtr cur_spec_ptr) {
  // get m/z and intensity values from the spectra
  std::vector<pwiz::msdata::MZIntensityPair> pairs;
  cur_spec_ptr->getMZIntensityPairs(pairs);
  LOG_DEBUG("mz pair size " << pairs.size());
  // make sure the peak list is already sorted
  std::sort(pairs.begin(), pairs.end(),
            [](const pwiz::msdata::MZIntensityPair &a, const pwiz::msdata::MZIntensityPair& b) {
            return a.mz < b.mz;
            });
  // get peaks
  PeakPtrVec peak_list; 
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i].intensity > 0.0) {
      PeakPtr peak_ptr = std::make_shared<Peak>(pairs[i].mz, pairs[i].intensity);
      peak_list.push_back(peak_ptr);
    }
  }
  return peak_list;
}

void PwMsReader::parsePrecursor(MsHeaderPtr header_ptr, 
                                pwiz::msdata::SpectrumInfo &spec_info,
                                pwiz::msdata::SpectrumPtr cur_spec_ptr) {
  // isolation window default values
  double prec_target_mz = 0;
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
  header_ptr->setPrecWinBegin(prec_target_mz - isolation_lower_offset);
  header_ptr->setPrecWinEnd(prec_target_mz + isolation_upper_offset);
  double prec_scan_num = -1;
  if (spec_info.precursors.size() > 0) {
    //get precursor scan ID from mzML
    prec_scan_num = parseNum(cur_spec_ptr->precursors[0].spectrumID, prev_ms1_scan_id_);
  }
  header_ptr->setMsOneScan(prec_scan_num);
}

void PwMsReader::parseActivation(MsHeaderPtr header_ptr, 
                                 pwiz::msdata::SpectrumInfo &spec_info,
                                 pwiz::msdata::SpectrumPtr cur_spec_ptr) {
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
    std::cout << "\nPlease select an activation method in database search.";
    std::cout << "\nExample: -a CID" << std::endl;
    exit(EXIT_FAILURE);
  }
  LOG_DEBUG("ac name " << ac_name);
  ActivationPtr activation_ptr = ActivationBase::getActivationPtrByName(ac_name);
  header_ptr->setActivationPtr(activation_ptr);
}

double PwMsReader::parseFaims(pwiz::msdata::SpectrumPtr cur_spec_ptr) {
  //add voltage information if it exists
  std::vector<pwiz::data::CVParam> cv_list = (*cur_spec_ptr).cvParams;
  for (size_t i = 0; i < cv_list.size(); i++) {
    if (cv_list[i].cvid == pwiz::cv::MS_FAIMS_compensation_voltage) {
      return std::stod(cv_list[i].value);
      break;
    }
  }
  return std::numeric_limits<double>::max();
}


bool PwMsReader::readOneMs(int sp_id, PeakPtrVec &peak_list, MsHeaderPtr &header_ptr) {
  bool get_binary_data = true;
  pwiz::msdata::SpectrumPtr cur_spec_ptr = spec_list_ptr_->spectrum(sp_id, get_binary_data);
  if (cur_spec_ptr == nullptr) {return false;}

  bool is_centroided = cur_spec_ptr->hasCVParam(pwiz::cv::MS_centroid_spectrum);
  if (!is_centroided) {
    std::cout << "Error: The data file contains profile data." << std::endl;
    std::cout << "TopFD can process only centroided, not profile, MS data." << std::endl;  
    exit(EXIT_FAILURE);
  }
  
  //get peaks;
  peak_list = parsePeaks(cur_spec_ptr);

  // scan info
  pwiz::msdata::SpectrumInfo spec_info(*cur_spec_ptr);

  header_ptr = std::make_shared<MsHeader>();
  header_ptr->setFileName(file_name_);
  parseScanNum(header_ptr, spec_info);
  header_ptr->setTitle("Scan_" + str_util::toString(spec_info.scanNumber));
  header_ptr->setRetentionTime(spec_info.retentionTime);
  double voltage = parseFaims(cur_spec_ptr);
  header_ptr->setVoltage(voltage); 

  int ms_level = spec_info.msLevel;
  LOG_DEBUG("ms_level " << ms_level);
  header_ptr->setMsLevel(ms_level);
  if (ms_level == 1) {
    //header_ptr->setSpecId(ms1_cnt_);
    prev_ms1_scan_id_ = spec_info.scanNumber;
    ms1_cnt_++;
  }
  else {
    //header_ptr->setSpecId(ms2_cnt_);
    parsePrecursor(header_ptr, spec_info, cur_spec_ptr);
    parseActivation(header_ptr, spec_info, cur_spec_ptr);
    ms2_cnt_++;
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

int PwMsReader::readNextWithVoltage(double voltage) {
  peak_list_.clear();
  header_ptr_ = nullptr;

  int found = false;
  while (!found) {
    if (input_sp_id_ >= input_sp_num_) {
      LOG_DEBUG("Only " << input_sp_num_  << " spectra in the input data!");
      return -1;
    }
    bool valid = readOneMs(input_sp_id_, peak_list_, header_ptr_); 
    if (valid && header_ptr_->getVoltage() == voltage) {
      found = true;
    }
    input_sp_id_++;
  }
  return 1;
}

MzmlProfilePtr PwMsReader::readProfile() {
  int ms_1_cnt = 0;
  int ms_2_cnt = 0;
  std::map<double, std::pair<int,int>> volt_map;
  bool get_binary_data = false;
  for (int sp_id = 0; sp_id < input_sp_num_; sp_id++) {
    pwiz::msdata::SpectrumPtr cur_spec_ptr = spec_list_ptr_->spectrum(sp_id, get_binary_data);
    if (cur_spec_ptr == nullptr) {continue;}
    pwiz::msdata::SpectrumInfo spec_info(*cur_spec_ptr);
    int ms_level = spec_info.msLevel;
    double voltage = parseFaims(cur_spec_ptr);
    if (ms_level == 1) {
      ms_1_cnt++;
    }
    else if (ms_level == 2) {
      ms_2_cnt++;
    }
    if (voltage != std::numeric_limits<double>::max()) {
      auto search = volt_map.find(voltage);
      // if not found
      if (search == volt_map.end()) {
        if (ms_level == 1) {
          std::pair<int,int> cnt(1,0);
          std::pair<double, std::pair<int,int>> new_volt(voltage, cnt);
          volt_map.insert(new_volt);
        }
        else if (ms_level == 2) {
          std::pair<int,int> cnt(0,1);
          std::pair<double, std::pair<int,int>> new_volt(voltage, cnt);
          volt_map.insert(new_volt);
        }
      }
      // if found
      else {
        if (ms_level == 1) {
          volt_map[voltage].first++;
        }
        else if (ms_level == 2) {
          volt_map[voltage].second++;
        }
      }
    }
  }
  MzmlProfilePtr profile_ptr = std::make_shared<MzmlProfile>(ms_1_cnt, ms_2_cnt, volt_map);
  return profile_ptr;
}

}  // namespace toppic
