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

#include <vector>
#include <string>
#include <algorithm>

#include "base/logger.hpp"
#include "base/activation_base.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/data/common/cv.hpp"
#include "pwiz/data/common/ParamTypes.hpp"
#include "feature/raw_ms_reader.hpp"

namespace prot {

RawMsReader::RawMsReader(const std::string & file_name):
    file_name_(file_name),
    input_sp_id_(0),
    output_sp_id_(0) {
      msd_ptr_ = std::make_shared<pwiz::msdata::MSDataFile>(file_name, &readers_);
      spec_list_ptr_ =  msd_ptr_->run.spectrumListPtr;
      input_sp_num_ = spec_list_ptr_->size();
    }

int RawMsReader::readNext() {
  peak_list_.clear();
  header_ptr_ = nullptr;

  if (input_sp_id_ >= input_sp_num_) {
    return -1;
  }

  pwiz::msdata::SpectrumPtr cur_spec_ptr = nullptr;
  // read m/z and intensity values from the spectra
  bool get_binary_data = true;
  while (cur_spec_ptr == nullptr) {
    cur_spec_ptr = spec_list_ptr_->spectrum(input_sp_id_, get_binary_data);
    input_sp_id_++;
    if (input_sp_id_ >= input_sp_num_ + 1) {
      LOG_ERROR("Only " << input_sp_num_  << " spectra in the input data!");
      return -1;
    }
  }

  std::vector<pwiz::msdata::MZIntensityPair> pairs;
  cur_spec_ptr->getMZIntensityPairs(pairs);
  LOG_DEBUG("mz pair size " << pairs.size());
  // make sure the peak list is already sorted
  std::sort(pairs.begin(), pairs.end(),
            [](const pwiz::msdata::MZIntensityPair &a, const pwiz::msdata::MZIntensityPair& b) {
            return a.mz < b.mz;
            });

  pwiz::msdata::SpectrumInfo spec_info(*cur_spec_ptr);
  peak_list_.empty();
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i].intensity > 0.0) {
      PeakPtr peak_ptr = std::make_shared<Peak>(pairs[i].mz, pairs[i].intensity);
      peak_list_.push_back(peak_ptr);
    }
  }
  int ms_level = spec_info.msLevel;
  LOG_DEBUG("ms_level " << ms_level);
  if (ms_level == 2) {
    double prec_mz;
    if (spec_info.precursors.size() == 0) {
      prec_mz = 0;
    } else {
      prec_mz = spec_info.precursors[0].mz;
    }

    if (prec_mz < 0) {
      prec_mz = 0;
    }

    int prec_charge;
    if (spec_info.precursors.size() == 0) {
      prec_charge = 1;
    } else {
      prec_charge = static_cast<int>(spec_info.precursors[0].charge);
    }

    if (prec_charge  < 0) {
      prec_charge = 1;
    }

    LOG_DEBUG("prec mz " << prec_mz << " scan number " << spec_info.scanNumber);
    header_ptr_ = std::make_shared<MsHeader>();
    header_ptr_->setId(ms2_cnt);
    ms2_cnt++;
    header_ptr_->setScan(spec_info.scanNumber);
    header_ptr_->setMsLevel(ms_level);
    header_ptr_->setPrecCharge(prec_charge);
    header_ptr_->setFileName(file_name_);
    header_ptr_->setTitle("Scan_" + std::to_string(spec_info.scanNumber));
    // here is average mz
    header_ptr_->setPrecSpMz(prec_mz);
    header_ptr_->setRetentionTime(spec_info.retentionTime);

    std::string ac_name;
    if (cur_spec_ptr->precursors.size() != 0) {
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
        }
      }
    }
    if (ac_name == "") {
      LOG_WARN("No activation information is available in reading the spectrum with scan " << spec_info.scanNumber);
    }
    LOG_DEBUG("ac name " << ac_name);
    ActivationPtr activation_ptr = ActivationBase::getActivationPtrByName(ac_name);
    header_ptr_->setActivationPtr(activation_ptr);
  } else {
    header_ptr_ = std::make_shared<MsHeader>();
    header_ptr_->setId(ms1_cnt);
    ms1_cnt++;
    header_ptr_->setScan(spec_info.scanNumber);
    header_ptr_->setMsLevel(ms_level);
    header_ptr_->setPrecCharge(0);
    header_ptr_->setFileName(file_name_);
    header_ptr_->setTitle("Scan_" + std::to_string(spec_info.scanNumber));
    header_ptr_->setRetentionTime(spec_info.retentionTime);
  }
  return 1;
}

}  // namespace prot
