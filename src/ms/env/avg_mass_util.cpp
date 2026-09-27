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

#include "ms/env/avg_mass_util.hpp"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/env/env_base.hpp"
#include "ms/spec/deconv_peak.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/spec/peak_util.hpp"
#include "ms/spec/precursor.hpp"

namespace toppic {

namespace avg_mass_util {

namespace {

// A precursor whose m/z is that of the given neutral mass (the other
// fields copied).
PrecursorPtr precursorWithMass(const PrecursorPtr& prec_ptr, double mass) {
  double mz = peak_util::compMz(mass, prec_ptr->getCharge());
  return std::make_shared<Precursor>(prec_ptr->getPrecId(),
                                     prec_ptr->getFeatureId(), mz,
                                     prec_ptr->getCharge(), prec_ptr->getInte());
}

// A copy of the spectrum with every mass converted by convert.
DeconvMsPtr convertMs(const DeconvMsPtr& ms_ptr, double (*convert)(double)) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  MsHeaderPtr new_header_ptr = std::make_shared<MsHeader>(*header_ptr);
  PrecursorPtrVec new_prec_ptr_vec;
  for (const PrecursorPtr& prec_ptr : header_ptr->getPrecPtrVec()) {
    double mass = convert(prec_ptr->getMonoMass());
    new_prec_ptr_vec.push_back(precursorWithMass(prec_ptr, mass));
  }
  new_header_ptr->setPrecPtrVec(new_prec_ptr_vec);
  DeconvPeakPtrVec new_peaks;
  for (size_t i = 0; i < ms_ptr->size(); i++) {
    DeconvPeakPtr peak_ptr = ms_ptr->getPeakPtr(i);
    double mass = convert(peak_ptr->getPosition());
    new_peaks.push_back(std::make_shared<DeconvPeak>(
        peak_ptr->getSpId(), peak_ptr->getPeakId(), mass,
        peak_ptr->getIntensity(), peak_ptr->getCharge(),
        peak_ptr->getScore()));
  }
  return std::make_shared<DeconvMs>(new_header_ptr, std::move(new_peaks));
}

double toAvgMass(double mono_mass) {
  // a mass beyond the distribution table is left unchanged
  if (EnvBase::getEnvByMonoMass(mono_mass) == nullptr) {
    return mono_mass;
  }
  return EnvBase::convertMonoMassToAvgMass(mono_mass);
}

double toMonoMass(double avg_mass) {
  return EnvBase::convertAvgMassToMonoMass(avg_mass);
}

}  // namespace

DeconvMsPtr getAvgMassMs(const DeconvMsPtr& ms_ptr) {
  return convertMs(ms_ptr, toAvgMass);
}

std::string monoMassFileName(const std::string& avg_sp_file_name) {
  std::string base_name = file_util::basename(avg_sp_file_name);
  if (str_util::endsWith(base_name, "_ms2")) {
    return base_name.substr(0, base_name.size() - 4) + "_mono_ms2.msalign";
  }
  return base_name + "_mono.msalign";
}

std::string convertToMonoMass(const std::string& avg_sp_file_name) {
  std::string mono_sp_file_name = monoMassFileName(avg_sp_file_name);
  MsAlignWriter writer(mono_sp_file_name);
  std::string para_str = msalign_util::readParaStr(avg_sp_file_name);
  if (!para_str.empty()) {
    para_str += "\n";
  }
  para_str +=
      "#Precursor and fragment masses converted from average masses to "
      "monoisotopic masses by TopPIC (--average-mass).";
  writer.writePara(para_str);
  MsAlignReader reader(avg_sp_file_name, 1, nullptr);
  DeconvMsPtr ms_ptr = reader.getNextMsPtr();
  int sp_num = 0;
  while (ms_ptr != nullptr) {
    writer.writeMs(convertMs(ms_ptr, toMonoMass));
    sp_num++;
    ms_ptr = reader.getNextMsPtr();
  }
  LOG_DEBUG(sp_num << " spectra converted to " << mono_sp_file_name);
  // the downstream steps look for the TopFD feature file under the new name
  std::string feature_file_name =
      file_util::basename(avg_sp_file_name) + ".feature";
  if (std::filesystem::exists(feature_file_name)) {
    bool overwrite = true;
    file_util::copyFile(feature_file_name,
                        file_util::basename(mono_sp_file_name) + ".feature",
                        overwrite);
  }
  return mono_sp_file_name;
}

}  // namespace avg_mass_util

}  // namespace toppic
