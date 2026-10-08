// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "ms/mzml/pw_ms_reader.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "common/base/activation_base.hpp"
#include "common/base/ion_type_base.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "pwiz/data/common/ParamTypes.hpp"
#include "pwiz/data/common/cv.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"

namespace toppic {

namespace {

const std::string& getPtcrName() {
  static const std::string name = "PTCR";
  return name;
}

// The activation of an MS/MS spectrum from the activation cvParams of its
// precursor, or "" when it has no known one. The set of terms decides, not
// their order: converters write the terms of a combined activation in
// different orders (msconvert: CID, ETD, ECD, PQD, HCD, then the
// supplemental term; ThermoRawFileParser: the first reaction, then the
// supplemental one), and neither writes the EThcD or ETciD term.
// - EThcD (MS:1002631), or ETD (MS:1000598) with beam-type CID (MS:1000422)
//   or supplemental beam-type CID (MS:1002678): ETHCD.
// - ETciD (MS:1003182), or ETD with CID (MS:1000133) or supplemental CID
//   (MS:1002679): ETCID. EThcD is tested first, so ETD with both kinds of
//   CID is ETHCD.
// - Otherwise ETD, HCD (beam-type CID, a kind of CID), CID, then MPD
//   (MS:1000435).
// The proton transfer terms are not dissociations (hasPtcrTerm).
std::string getActivationName(const std::vector<pwiz::data::CVParam>& cv_list) {
  bool cid = false;
  bool hcd = false;
  bool etd = false;
  bool mpd = false;
  bool supp_cid = false;
  bool supp_hcd = false;
  bool ethcd = false;
  bool etcid = false;
  for (size_t i = 0; i < cv_list.size(); i++) {
    LOG_DEBUG("cv list " << i << " " << cv_list[i].cvid << " "
                         << cv_list[i].value);
    pwiz::cv::CVID cvid = cv_list[i].cvid;
    cid = cid || cvid == pwiz::cv::MS_CID;
    hcd = hcd || cvid == pwiz::cv::MS_HCD;
    etd = etd || cvid == pwiz::cv::MS_ETD;
    mpd = mpd || cvid == pwiz::cv::MS_MPD;
    supp_cid = supp_cid ||
               cvid == pwiz::cv::MS_supplemental_collision_induced_dissociation;
    supp_hcd =
        supp_hcd ||
        cvid ==
            pwiz::cv::MS_supplemental_beam_type_collision_induced_dissociation;
    ethcd = ethcd || cvid == pwiz::cv::MS_EThcD;
    etcid = etcid || cvid == pwiz::cv::MS_ETciD;
  }
  if (ethcd || (etd && (hcd || supp_hcd))) {
    return "ETHCD";
  }
  if (etcid || (etd && (cid || supp_cid))) {
    return "ETCID";
  }
  if (etd) {
    return "ETD";
  }
  if (hcd) {
    return "HCD";
  }
  if (cid) {
    return "CID";
  }
  if (mpd) {
    return "MPD";
  }
  return "";
}

// Whether activation terms have proton transfer charge reduction
// (MS:1003249) or proton transfer reaction (MS:1003248)
bool hasPtcrTerm(const std::vector<pwiz::data::CVParam>& cv_list) {
  for (const pwiz::data::CVParam& cv : cv_list) {
    if (cv.cvid == pwiz::cv::MS_proton_transfer_charge_reduction ||
        cv.cvid == pwiz::cv::MS_proton_transfer_reaction) {
      return true;
    }
  }
  return false;
}

// One stage of an MSn scan in its Thermo filter string: the isolation target
// and the activation terms
struct FilterStage {
  double target = 0;
  std::vector<pwiz::data::CVParam> terms;
};

// The stages of an MS/MS scan in its Thermo filter string (MS:1000512), in
// order, each with its target and the terms of its activations (cid, hcd,
// etd, ptr; others have none); none when the filter string does not list
// ms_level - 1 stages.
std::vector<FilterStage> getFilterStages(
    const pwiz::msdata::SpectrumPtr& spec, int ms_level) {
  std::vector<FilterStage> stages;
  if (ms_level < 2 || spec->scanList.scans.empty()) {
    return stages;
  }
  std::string filter =
      spec->scanList.scans[0].cvParam(pwiz::cv::MS_filter_string).value;
  std::string key = "ms" + std::to_string(ms_level) + " ";
  std::size_t pos = filter.find(key);
  std::size_t end =
      pos == std::string::npos ? std::string::npos : filter.find('[', pos);
  if (end == std::string::npos) {
    return stages;
  }
  std::stringstream ss(
      filter.substr(pos + key.size(), end - pos - key.size()));
  std::string item;
  while (ss >> item) {
    // <target>@<type><energy>[@<type><energy>...]
    std::vector<std::string> parts = str_util::split(item, "@");
    FilterStage stage;
    try {
      stage.target = std::stod(parts[0]);
    } catch (const std::exception&) {
      return std::vector<FilterStage>();
    }
    for (std::size_t i = 1; i < parts.size(); i++) {
      std::string type;
      for (char c : parts[i]) {
        if (!std::isalpha(static_cast<unsigned char>(c))) {
          break;
        }
        type += c;
      }
      if (type == "cid") {
        stage.terms.emplace_back(pwiz::cv::MS_CID);
      } else if (type == "hcd") {
        stage.terms.emplace_back(pwiz::cv::MS_HCD);
      } else if (type == "etd") {
        stage.terms.emplace_back(pwiz::cv::MS_ETD);
      } else if (type == "ptr") {
        stage.terms.emplace_back(
            pwiz::cv::MS_proton_transfer_charge_reduction);
      }
    }
    stages.push_back(stage);
  }
  if (static_cast<int>(stages.size()) != ms_level - 1) {
    stages.clear();
  }
  return stages;
}

// The stages of an MS/MS scan whose converter wrote fewer precursors than
// stages (getFilterStages), none otherwise. ThermoRawFileParser 2.0 writes
// one precursor for the MS3 scan "FTMS + p ESI Full ms3 1402.5000@cid60.00
// 1100.0000@ptr10.00 [...]", the last stage (its isolation window, the
// proton transfer term and the scan's charge state), and loses the
// dissociation of the first stage; and none for a scan whose precursor scan
// is not in the file (a scan cut out of a longer run), e.g. the MS2 scan
// "... Full ms2 1060.1759@ptr5.00 [...]".
std::vector<FilterStage> getMissingStages(
    const pwiz::msdata::SpectrumPtr& spec, int ms_level) {
  if (static_cast<int>(spec->precursors.size()) >= ms_level - 1) {
    return std::vector<FilterStage>();
  }
  return getFilterStages(spec, ms_level);
}

// The isolation window of a precursor: target, lower and upper offset (-1
// each when the precursor has none)
void getIsolationWindow(const pwiz::msdata::Precursor& prec, double& target,
                        double& lower, double& upper) {
  target = -1;
  lower = -1;
  upper = -1;
  for (const pwiz::data::CVParam& cv : prec.isolationWindow.cvParams) {
    if (cv.cvid == pwiz::cv::MS_isolation_window_target_m_z) {
      target = std::stod(cv.value);
    } else if (cv.cvid == pwiz::cv::MS_isolation_window_lower_offset) {
      lower = std::stod(cv.value);
    } else if (cv.cvid == pwiz::cv::MS_isolation_window_upper_offset) {
      upper = std::stod(cv.value);
    }
  }
}

// One stage of an MS/MS scan in acquisition order: its isolation target and
// window (-1 when unknown) and its activation terms
struct ScanStage {
  double target = -1;
  double lower = -1;
  double upper = -1;
  std::vector<pwiz::data::CVParam> terms;
};

// The stages of an MS/MS scan in acquisition order: the stages of the filter
// string when it lists ms_level - 1 of them (getFilterStages), each with the
// isolation window and activation terms of the file's precursor of the same
// target too (msconvert writes one precursor per stage, ThermoRawFileParser
// one for the last stage of an MS3 scan); else the file's precursors in file
// order (msconvert writes the first stage first)
std::vector<ScanStage> getScanStages(const pwiz::msdata::SpectrumPtr& spec,
                                     int ms_level) {
  std::vector<ScanStage> stages;
  std::vector<FilterStage> filter_stages = getFilterStages(spec, ms_level);
  if (filter_stages.empty()) {
    for (const pwiz::msdata::Precursor& prec : spec->precursors) {
      ScanStage stage;
      getIsolationWindow(prec, stage.target, stage.lower, stage.upper);
      stage.terms = prec.activation.cvParams;
      stages.push_back(stage);
    }
    return stages;
  }
  for (const FilterStage& filter_stage : filter_stages) {
    ScanStage stage;
    stage.target = filter_stage.target;
    stage.terms = filter_stage.terms;
    for (const pwiz::msdata::Precursor& prec : spec->precursors) {
      double target = 0;
      double lower = 0;
      double upper = 0;
      getIsolationWindow(prec, target, lower, upper);
      if (target > 0 && std::abs(target - filter_stage.target) <= 0.01) {
        stage.lower = lower;
        stage.upper = upper;
        stage.terms.insert(stage.terms.end(), prec.activation.cvParams.begin(),
                           prec.activation.cvParams.end());
        break;
      }
    }
    stages.push_back(stage);
  }
  return stages;
}

// The dissociation stage of an MS/MS scan that comes after its last proton
// transfer stage (the dissociation of a charge-reduced precursor, e.g. the
// MS3 scan "X@ptr5 Y@hcd80"), -1 when there is none: no proton transfer, no
// dissociation, or the proton transfer after the last dissociation or in its
// stage (a fragment PTCR spectrum, "X@cid60 Y@ptr10")
int getReducedPrecStage(const std::vector<ScanStage>& stages) {
  int last_ptcr = -1;
  int last_diss = -1;
  for (std::size_t i = 0; i < stages.size(); i++) {
    if (hasPtcrTerm(stages[i].terms)) {
      last_ptcr = static_cast<int>(i);
    }
    if (getActivationName(stages[i].terms) != "") {
      last_diss = static_cast<int>(i);
    }
  }
  if (last_ptcr < 0 || last_diss <= last_ptcr) {
    return -1;
  }
  return last_diss;
}

// The activation of PTCR spectra. It is not in activation_base.xml: a PTCR
// spectrum has no fragments, TopFD writes it to its own msalign file and no
// search reads it. The ion types are placeholders.
ActivationPtr getPtcrActivationPtr() {
  static const ActivationPtr ptcr_ptr = std::make_shared<Activation>(
      getPtcrName(), IonTypeBase::getIonTypePtrByName("B"),
      IonTypeBase::getIonTypePtrByName("Y"));
  return ptcr_ptr;
}

}  // namespace

PwMsReader::PwMsReader(const std::string& file_name) { init(file_name); }

PwMsReader::PwMsReader(const std::string& file_name, double isolation_window) {
  init(file_name);
  isolation_window_ = isolation_window;
}

PwMsReader::PwMsReader(const std::string& file_name, double isolation_window,
                       const std::string& activation) {
  init(file_name);
  isolation_window_ = isolation_window;
  activation_ = activation;
}

void PwMsReader::init(const std::string& file_name) {
  file_name_ = file_name;
  input_sp_id_ = 0;
  try {
    msd_ptr_ =
        std::make_shared<pwiz::msdata::MSDataFile>(file_name, &readers_);
  } catch (const std::exception& e) {
    // pwiz reports a truncated mzML file (missing index and closing tags)
    // as a bad input stream after its fallback re-indexing parse hits the
    // end of the file.
    LOG_ERROR("Unable to read the spectrum file " << file_name << ": "
              << e.what());
    std::cout << "\nERROR: Unable to read the spectrum file " << file_name
              << ".";
    std::string what = e.what();
    if (what.find("Bad istream") != std::string::npos) {
      std::cout << "\nThe mzML file is truncated or incomplete: it ends before "
                << "its closing tags and index.";
      std::cout << "\nCopy or convert the file again and check that it ends "
                << "with </indexedmzML> or </mzML>.";
    } else {
      std::cout << "\nReader message: " << what;
    }
    std::cout << std::endl;
    exit(EXIT_FAILURE);
  }
  spec_list_ptr_ = msd_ptr_->run.spectrumListPtr;
  input_sp_num_ = spec_list_ptr_->size();
  is_waters_instrument_ = checkWatersInstrument();
}

bool PwMsReader::checkWatersInstrument() {
  for (size_t i = 0; i < msd_ptr_->instrumentConfigurationPtrs.size(); i++) {
    pwiz::msdata::InstrumentConfigurationPtr conf_ptr =
        msd_ptr_->instrumentConfigurationPtrs[i];
    pwiz::msdata::CVParam param =
        conf_ptr->cvParam(pwiz::msdata::CVID::MS_Waters_instrument_model);
    if (!param.empty()) {
      LOG_DEBUG("Found waters model");
      return true;
    }
  }
  return false;
}

bool PwMsReader::checkCentroidData() {
  std::vector<pwiz::msdata::DataProcessingPtr> proc_ptr_vec =
      msd_ptr_->dataProcessingPtrs;
  for (size_t i = 0; i < proc_ptr_vec.size(); i++) {
    for (size_t j = 0; j < proc_ptr_vec[i]->processingMethods.size(); j++) {
      pwiz::msdata::ProcessingMethod proc_method =
          proc_ptr_vec[i]->processingMethods[j];
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

int PwMsReader::parseNum(const std::string& id, int default_scan) {
  std::string delimiter = "=";
  //
  if (id.substr(id.find_last_of(delimiter) + 1) == "") {
    LOG_INFO("Scan number information is missing!");
    return default_scan;
  } else {
    return std::stoi(id.substr(id.find_last_of(delimiter) + 1));
  }
}

void PwMsReader::parseScanNum(const MsHeaderPtr& header_ptr,
                              pwiz::msdata::SpectrumInfo& spec_info) {
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
            [](const pwiz::msdata::MZIntensityPair& a,
               const pwiz::msdata::MZIntensityPair& b) { return a.mz < b.mz; });
  // get peaks
  PeakPtrVec peak_list;
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i].intensity > 0.0 && pairs[i].intensity < MAX_INTE_ &&
        pairs[i].mz > 0.0 && pairs[i].mz < MAX_MZ_) {
      PeakPtr peak_ptr =
          std::make_shared<Peak>(pairs[i].mz, pairs[i].intensity);
      peak_list.push_back(peak_ptr);
    }
  }
  return peak_list;
}

void PwMsReader::parsePrecursor(const MsHeaderPtr& header_ptr,
                                pwiz::msdata::SpectrumInfo& spec_info,
                                pwiz::msdata::SpectrumPtr cur_spec_ptr) {
  // isolation window default values
  double prec_target_mz = 0;
  double isolation_lower_offset = isolation_window_ / 2;
  double isolation_upper_offset = isolation_window_ / 2;
  if (cur_spec_ptr->precursors.size() > 0) {
    std::vector<pwiz::data::CVParam> cv_list =
        cur_spec_ptr->precursors[0].isolationWindow.cvParams;
    for (size_t i = 0; i < cv_list.size(); i++) {
      LOG_DEBUG("cv list " << i << " " << cv_list[i].cvid << " "
                           << cv_list[i].value);
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
  // a scan with fewer precursors than stages (getMissingStages): the file's
  // precursor, if any, is a later stage when its target is not the first
  // stage's, and the precursor window is then the first stage's target with
  // the default window (-w)
  std::vector<FilterStage> stages =
      getMissingStages(cur_spec_ptr, spec_info.msLevel);
  bool first_stage_missing =
      !stages.empty() && std::abs(prec_target_mz - stages[0].target) > 0.01;
  if (first_stage_missing) {
    prec_target_mz = stages[0].target;
    isolation_lower_offset = isolation_window_ / 2;
    isolation_upper_offset = isolation_window_ / 2;
  }
  LOG_DEBUG("Precursor target mz " << prec_target_mz << " lower offset "
                                   << isolation_lower_offset << " upper offset "
                                   << isolation_upper_offset);
  header_ptr->setPrecTargetMz(prec_target_mz);
  header_ptr->setPrecWinBegin(prec_target_mz - isolation_lower_offset);
  header_ptr->setPrecWinEnd(prec_target_mz + isolation_upper_offset);
  // the selected ion as the file gives it (used for PTCR spectra)
  if (cur_spec_ptr->precursors.size() > 0 &&
      !cur_spec_ptr->precursors[0].selectedIons.empty()) {
    const pwiz::msdata::SelectedIon& ion =
        cur_spec_ptr->precursors[0].selectedIons[0];
    for (const pwiz::data::CVParam& cv : ion.cvParams) {
      if (cv.value.empty()) {
        continue;
      }
      if (cv.cvid == pwiz::cv::MS_selected_ion_m_z) {
        header_ptr->setSelectedMz(std::stod(cv.value));
      } else if (cv.cvid == pwiz::cv::MS_charge_state) {
        header_ptr->setSelectedCharge(std::stoi(cv.value));
      }
    }
  }
  if (first_stage_missing) {
    // the first stage's target; the charge state stays the file's (the
    // scan's, ThermoRawFileParser)
    header_ptr->setSelectedMz(stages[0].target);
  }
  // the isolation window of the proton transfer stage: that of the last
  // precursor with a proton transfer term (the second precursor of an MS3
  // scan "X@cid Y@ptr"), else the precursor window
  header_ptr->setPtcrWin(header_ptr->getPrecWinBegin(),
                         header_ptr->getPrecWinEnd());
  // missing stages: the target of the last proton transfer stage of the
  // filter string with the default window, unless a precursor of the file
  // is that stage (below)
  for (const FilterStage& stage : stages) {
    if (hasPtcrTerm(stage.terms)) {
      header_ptr->setPtcrWin(stage.target - isolation_window_ / 2,
                             stage.target + isolation_window_ / 2);
    }
  }
  for (const pwiz::msdata::Precursor& prec : cur_spec_ptr->precursors) {
    if (!hasPtcrTerm(prec.activation.cvParams)) {
      continue;
    }
    double target = 0;
    double lower = isolation_window_ / 2;
    double upper = isolation_window_ / 2;
    for (const pwiz::data::CVParam& cv : prec.isolationWindow.cvParams) {
      if (cv.cvid == pwiz::cv::MS_isolation_window_target_m_z) {
        target = std::stod(cv.value);
      } else if (cv.cvid == pwiz::cv::MS_isolation_window_lower_offset) {
        lower = std::stod(cv.value);
      } else if (cv.cvid == pwiz::cv::MS_isolation_window_upper_offset) {
        upper = std::stod(cv.value);
      }
    }
    if (target > 0) {
      header_ptr->setPtcrWin(target - lower, target + upper);
    }
  }
  double prec_scan_num = -1;
  if (spec_info.precursors.size() > 0) {
    // get precursor scan ID from mzML
    prec_scan_num =
        parseNum(cur_spec_ptr->precursors[0].spectrumID, prev_ms1_scan_id_);
  }
  if (spec_info.msLevel > 2) {
    // an MSn scan belongs to the MS1 scan of its cycle; its precursor
    // reference is the MS(n-1) scan
    prec_scan_num = prev_ms1_scan_id_;
    // the MS(n-1) scan: the spectrumRef of the precursor of the last stage
    // of the filter string (ThermoRawFileParser writes the last stage
    // first, msconvert the first stage first), else of the last precursor
    if (!cur_spec_ptr->precursors.empty()) {
      const pwiz::msdata::Precursor* parent = &cur_spec_ptr->precursors.back();
      std::vector<FilterStage> filter_stages =
          getFilterStages(cur_spec_ptr, spec_info.msLevel);
      if (!filter_stages.empty()) {
        for (const pwiz::msdata::Precursor& prec : cur_spec_ptr->precursors) {
          double target = 0;
          double lower = 0;
          double upper = 0;
          getIsolationWindow(prec, target, lower, upper);
          if (target > 0 &&
              std::abs(target - filter_stages.back().target) <= 0.01) {
            parent = &prec;
            break;
          }
        }
      }
      header_ptr->setParentScan(parseNum(parent->spectrumID, -1));
    }
  }
  header_ptr->setMsOneScan(prec_scan_num);
  // the dissociation stage of a charge-reduced precursor (an MSn scan "X@ptr
  // Y@hcd"): its target and window, the default window (-w) when the file
  // has none
  std::vector<ScanStage> scan_stages =
      getScanStages(cur_spec_ptr, spec_info.msLevel);
  int reduced_stage = getReducedPrecStage(scan_stages);
  if (reduced_stage >= 0) {
    const ScanStage& stage = scan_stages[reduced_stage];
    double lower = stage.lower >= 0 ? stage.lower : isolation_window_ / 2;
    double upper = stage.upper >= 0 ? stage.upper : isolation_window_ / 2;
    header_ptr->setReducedPrecWin(stage.target, stage.target - lower,
                                  stage.target + upper);
  }
}

void PwMsReader::parseActivation(const MsHeaderPtr& header_ptr,
                                 pwiz::msdata::SpectrumInfo& spec_info,
                                 pwiz::msdata::SpectrumPtr cur_spec_ptr) {
  // the activation terms of all precursors: an MSn scan has one precursor
  // per stage (msconvert writes the MS2 stage first), e.g. CID of the
  // precursor, then proton transfer charge reduction of its fragments
  std::vector<pwiz::data::CVParam> terms;
  for (const pwiz::msdata::Precursor& prec : cur_spec_ptr->precursors) {
    terms.insert(terms.end(), prec.activation.cvParams.begin(),
                 prec.activation.cvParams.end());
  }
  // and those of the stages the converter did not write (getMissingStages)
  for (const FilterStage& stage :
       getMissingStages(cur_spec_ptr, spec_info.msLevel)) {
    terms.insert(terms.end(), stage.terms.begin(), stage.terms.end());
  }
  bool ptcr = hasPtcrTerm(terms);
  std::string ac_name = activation_;
  if (ac_name == getPtcrName()) {
    // topfd --ptcr: every MS/MS spectrum is a PTCR spectrum
    header_ptr->setPtcr(true);
    header_ptr->setActivationPtr(getPtcrActivationPtr());
    return;
  }
  if (ac_name == "" || ac_name == "FILE") {
    ac_name = getActivationName(terms);
    if (ac_name == "" && ptcr) {
      // proton transfer without dissociation: a PTCR spectrum of the
      // isolated ions
      header_ptr->setPtcr(true);
      header_ptr->setActivationPtr(getPtcrActivationPtr());
      return;
    }
  }
  if (ac_name == "") {
    LOG_WARN(
        "No activation information is available in reading the spectrum with "
        "scan "
        << spec_info.scanNumber);
    std::cout << "\nERROR: Unable to read the activation method from the file.";
    std::cout << "\nPlease select an activation method in database search.";
    std::cout << "\nExample: -a CID" << std::endl;
    exit(EXIT_FAILURE);
  }
  LOG_DEBUG("ac name " << ac_name);
  ActivationPtr activation_ptr =
      ActivationBase::getActivationPtrByName(ac_name);
  header_ptr->setActivationPtr(activation_ptr);
  // a dissociation (from the file or -a) and a proton transfer term: the
  // fragments, charge-reduced (a fragment PTCR spectrum), unless the
  // dissociation stage comes after the proton transfer stage: then the
  // dissociation of a charge-reduced precursor (getReducedPrecStage)
  bool reduced_prec =
      ptcr && getReducedPrecStage(getScanStages(cur_spec_ptr,
                                                spec_info.msLevel)) >= 0;
  header_ptr->setReducedPrec(reduced_prec);
  header_ptr->setFragPtcr(ptcr && !reduced_prec);
}

double PwMsReader::parseFaims(pwiz::msdata::SpectrumPtr cur_spec_ptr) {
  // add voltage information if it exists
  std::vector<pwiz::data::CVParam> cv_list = (*cur_spec_ptr).cvParams;
  for (size_t i = 0; i < cv_list.size(); i++) {
    if (cv_list[i].cvid == pwiz::cv::MS_FAIMS_compensation_voltage) {
      return std::stod(cv_list[i].value);
      break;
    }
  }
  return std::numeric_limits<double>::max();
}

bool PwMsReader::readOneMs(int sp_id, PeakPtrVec& peak_list,
                           MsHeaderPtr& header_ptr) {
  bool get_binary_data = true;
  pwiz::msdata::SpectrumPtr cur_spec_ptr =
      spec_list_ptr_->spectrum(sp_id, get_binary_data);
  if (cur_spec_ptr == nullptr) {
    return false;
  }

  bool is_centroided = cur_spec_ptr->hasCVParam(pwiz::cv::MS_centroid_spectrum);
  if (!is_centroided) {
    std::cout << "Error: The data file contains profile data." << std::endl;
    std::cout << "TopFD can process only centroided, not profile, MS data."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // get peaks;
  peak_list = parsePeaks(cur_spec_ptr);

  // scan info
  pwiz::msdata::SpectrumInfo spec_info(*cur_spec_ptr);

  header_ptr = std::make_shared<MsHeader>();
  header_ptr->setFileName(file_name_);
  parseScanNum(header_ptr, spec_info);
  header_ptr->setTitle("Scan_" + std::to_string(spec_info.scanNumber));
  header_ptr->setRetentionTime(spec_info.retentionTime);
  double voltage = parseFaims(cur_spec_ptr);
  header_ptr->setVoltage(voltage);

  int ms_level = spec_info.msLevel;
  LOG_DEBUG("ms_level " << ms_level);
  header_ptr->setMsLevel(ms_level);
  if (ms_level == 1) {
    // header_ptr->setSpecId(ms1_cnt_);
    prev_ms1_scan_id_ = spec_info.scanNumber;
    ms1_cnt_++;
  } else {
    // header_ptr->setSpecId(ms2_cnt_);
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
      LOG_DEBUG("Only " << input_sp_num_ << " spectra in the input data!");
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
      LOG_DEBUG("Only " << input_sp_num_ << " spectra in the input data!");
      // The last spectrum read while skipping other voltages is still in
      // header_ptr_/peak_list_; drop it, or the caller (which only tests
      // header_ptr_ for null) would emit it into this voltage's outputs.
      peak_list_.clear();
      header_ptr_ = nullptr;
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

namespace {

// Returns true if the spectrum's precursor specifies its isolation window width
// (a lower/upper isolation-window offset cvParam) -- i.e. the file carries the
// MS/MS precursor window rather than relying on a default width.
bool specHasPrecWindow(const pwiz::msdata::SpectrumPtr& spec) {
  if (spec->precursors.empty()) {
    return false;
  }
  for (const pwiz::data::CVParam& cv :
       spec->precursors[0].isolationWindow.cvParams) {
    if (cv.cvid == pwiz::cv::MS_isolation_window_lower_offset ||
        cv.cvid == pwiz::cv::MS_isolation_window_upper_offset) {
      return true;
    }
  }
  return false;
}

// Whether a spectrum is an MS/MS spectrum TopFD reads
// (MzmlMsGroupReader): an MS2 scan, or an MSn scan (n > 2) with a proton
// transfer term, of a precursor or of a stage the converter did not write
// (all of them with topfd --ptcr)
bool isReadMsTwo(const pwiz::msdata::SpectrumPtr& spec, int ms_level,
                 bool ptcr_forced) {
  if (ms_level == 2) {
    return true;
  }
  if (ms_level < 2) {
    return false;
  }
  if (ptcr_forced) {
    return true;
  }
  for (const pwiz::msdata::Precursor& prec : spec->precursors) {
    if (hasPtcrTerm(prec.activation.cvParams)) {
      return true;
    }
  }
  for (const FilterStage& stage : getMissingStages(spec, ms_level)) {
    if (hasPtcrTerm(stage.terms)) {
      return true;
    }
  }
  return false;
}

}  // namespace

MzmlProfilePtr PwMsReader::readProfile() {
  int ms_1_cnt = 0;
  int ms_2_cnt = 0;
  int ms_3_cnt = 0;
  bool all_ms2_have_prec_window = true;
  std::map<double, std::pair<int, int>> volt_map;
  bool get_binary_data = false;
  for (int sp_id = 0; sp_id < input_sp_num_; sp_id++) {
    pwiz::msdata::SpectrumPtr cur_spec_ptr =
        spec_list_ptr_->spectrum(sp_id, get_binary_data);
    if (cur_spec_ptr == nullptr) {
      continue;
    }
    pwiz::msdata::SpectrumInfo spec_info(*cur_spec_ptr);
    int ms_level = spec_info.msLevel;
    // the MS/MS spectra TopFD reads (MS3 scans with a proton transfer stage
    // too); other MSn scans are not counted
    if (ms_level > 2 &&
        !isReadMsTwo(cur_spec_ptr, ms_level, activation_ == getPtcrName())) {
      // an MS3 scan without a proton transfer stage: DeconvMs3Process
      if (ms_level == 3) ms_3_cnt++;
      continue;
    }
    double voltage = parseFaims(cur_spec_ptr);
    if (ms_level == 1) {
      ms_1_cnt++;
    } else if (ms_level >= 2) {
      ms_2_cnt++;
      if (all_ms2_have_prec_window && !specHasPrecWindow(cur_spec_ptr)) {
        all_ms2_have_prec_window = false;
      }
    }
    if (voltage != std::numeric_limits<double>::max()) {
      auto search = volt_map.find(voltage);
      // if not found
      if (search == volt_map.end()) {
        if (ms_level == 1) {
          std::pair<int, int> cnt(1, 0);
          std::pair<double, std::pair<int, int>> new_volt(voltage, cnt);
          volt_map.insert(new_volt);
        } else if (ms_level >= 2) {
          std::pair<int, int> cnt(0, 1);
          std::pair<double, std::pair<int, int>> new_volt(voltage, cnt);
          volt_map.insert(new_volt);
        }
      }
      // if found
      else {
        if (ms_level == 1) {
          volt_map[voltage].first++;
        } else if (ms_level >= 2) {
          volt_map[voltage].second++;
        }
      }
    }
  }
  // The file carries the MS/MS precursor windows only when every MS/MS scan
  // specifies one (and there is at least one MS/MS scan).
  bool has_prec_window = (ms_2_cnt > 0) && all_ms2_have_prec_window;
  MzmlProfilePtr profile_ptr = std::make_shared<MzmlProfile>(
      ms_1_cnt, ms_2_cnt, has_prec_window, volt_map);
  profile_ptr->setMs3Cnt(ms_3_cnt);
  return profile_ptr;
}

}  // namespace toppic
