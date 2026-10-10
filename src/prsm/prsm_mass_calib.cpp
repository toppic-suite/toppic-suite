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

#include "prsm/prsm_mass_calib.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/deconv_peak.hpp"
#include "ms/spec/ms_header.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/theo_peak_util.hpp"
#include "seq/alter_type.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

namespace prsm_mass_calib {

MassCalibrantSet readCalibrants(const PrsmParaPtr& prsm_para_ptr,
                                const std::string& input_file_ext,
                                double max_evalue, double tol_ppm) {
  MassCalibrantSet set;
  PrsmPtrVec prsm_ptrs =
      prsm_reader_util::readPrsmsWithSpectra(prsm_para_ptr, input_file_ext);
  double min_mass = prsm_para_ptr->getSpParaPtr()->getMinMass();
  set.prsm_num = static_cast<int>(prsm_ptrs.size());
  for (const PrsmPtr& prsm_ptr : prsm_ptrs) {
    ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
    if (form_ptr->getSeqName().find("DECOY_") == 0 ||
        !(prsm_ptr->getEValue() <= max_evalue) ||
        form_ptr->getAlterNum(AlterType::UNEXPECTED) != 0 ||
        form_ptr->getVarPtmNum() != 0) {
      continue;
    }
    set.used_prsm_num++;
    DeconvMsPtr ms_ptr = prsm_ptr->getDeconvMsPtrVec()[0];
    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    int scan = header_ptr->getFirstScanNum();
    // MsHeader keeps the retention time in seconds
    double rt_min = header_ptr->getRetentionTime() / 60;
    mass_calib_fit::addPrecCalibrant(scan, rt_min,
                                     header_ptr->getFirstPrecCharge(),
                                     prsm_ptr->getOriPrecMass(),
                                     form_ptr->getMass(), &set);
    TheoPeakPtrVec theo_peaks = theo_peak_util::geneProteoformTheoPeak(
        form_ptr, header_ptr->getActivationPtr(), min_mass,
        prsm_ptr->getNTermLabelMass());
    std::vector<double> theo_masses = theo_peak_util::getTheoMassVec(theo_peaks);
    std::sort(theo_masses.begin(), theo_masses.end());
    std::vector<std::pair<double, int>> masses;
    for (const DeconvPeakPtr& peak_ptr : ms_ptr->getPeakPtrVec()) {
      masses.emplace_back(peak_ptr->getMonoMass(), peak_ptr->getCharge());
    }
    mass_calib_fit::addFragCalibrants(scan, rt_min, theo_masses, masses,
                                      tol_ppm, &set);
  }
  return set;
}

namespace {

// <dir>/<base> for <dir>/<base>_ms2.msalign
std::string fileBase(const std::string& sp_file_name) {
  std::string base = file_util::basename(sp_file_name);
  if (str_util::endsWith(base, "_ms2")) {
    base = base.substr(0, base.size() - 4);
  }
  return base;
}

}  // namespace

std::string calibSpFileName(const std::string& sp_file_name) {
  return fileBase(sp_file_name) + "_calib_ms2.msalign";
}

std::string modelFileName(const std::string& sp_file_name) {
  return fileBase(sp_file_name) + "_mass_calib.json";
}

}  // namespace prsm_mass_calib

}  // namespace toppic
