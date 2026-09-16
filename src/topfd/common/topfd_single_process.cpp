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
//
#include "topfd/common/topfd_single_process.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "common/base/activation_base.hpp"
#include "common/base/base_data.hpp"
#include "common/util/str_util.hpp"
#include "common/util/time_util.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/env/match_env_writer.hpp"
#include "ms/mzml/mzml_ms.hpp"
#include "ms/mzml/mzml_ms_sql_writer.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "topfd/deconv/deconv_single_sp.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

namespace topfd_single_process {

namespace {

PeakPtrVec readPeakFile(const std::string& file_name) {
  PeakPtrVec peak_list;
  std::ifstream input;
  input.open(file_name.c_str(), std::ios::in);
  std::string line;
  while (std::getline(input, line)) {
    str_util::trim(line);
    if (line.length() == 0) {
      continue;
    }
    std::vector<std::string> strs = str_util::split(line, " ");
    double mz = std::stod(strs[0]);
    double inte = std::stod(strs[1]);
    PeakPtr peak_ptr = std::make_shared<Peak>(mz, inte);
    peak_list.push_back(peak_ptr);
  }
  input.close();
  // Deconvolution assumes peaks in increasing m/z order; the input file need
  // not be sorted. A stable sort keeps the file order of equal m/z values.
  std::stable_sort(peak_list.begin(), peak_list.end(), Peak::cmpPosInc);
  return peak_list;
}

void processOneFile(const TopfdParaPtr& para_ptr,
                    const std::string& spec_file_name) {
  try {
    int ms_level = 2;
    double max_mass = para_ptr->getMaxMass();
    double max_charge = para_ptr->getMaxCharge();

    // The text-peak-list input is a single MS/MS spectrum with no MS1 scan.
    // createSqlDb (called below through setMzmlFileNameAndFaims) writes these
    // into the ms_info table, so set them first (they default to -1, which is
    // only updated in the mzML flow).
    para_ptr->setMs1ScanNumber(0);
    para_ptr->setMs2ScanNumber(1);
    // Name the outputs after the input file, as the mzML flow does: the base
    // name is the input path minus its extension, so <input>_ms2.msalign,
    // <input>_ms2.env and <input>.sqlite land next to the input file. This
    // also creates the SQLite database when it is enabled.
    para_ptr->setMzmlFileNameAndFaims(spec_file_name, false, -1);

    PeakPtrVec peak_list = readPeakFile(spec_file_name);

    MatchEnvPtrVec result_envs;
    if (peak_list.size() > 0) {
      DeconvSingleSpPtr deconv_ptr = std::make_shared<DeconvSingleSp>(
          para_ptr, peak_list, ms_level, max_mass, max_charge);
      result_envs = deconv_ptr->deconv();
    }

    // header
    MsHeaderPtr header_ptr = std::make_shared<MsHeader>();
    header_ptr->setSpecId(0);
    header_ptr->setSingleScan(1);

    DeconvMsPtr ms_ptr =
        match_env_util::getDeconvMsPtr(header_ptr, result_envs);

    std::string output_base_name = para_ptr->getOutputBaseName();
    std::string ms2_msalign_name = output_base_name + "_ms2.msalign";
    MsAlignWriterPtr ms2_writer_ptr =
        std::make_shared<MsAlignWriter>(ms2_msalign_name);
    ms2_writer_ptr->writeMs(ms_ptr);
    ms2_writer_ptr = nullptr;
    std::string ms_env_name = output_base_name + "_ms2.env";
    match_env_writer::writePeakList(ms_env_name, peak_list, result_envs);

    // Optionally write the deconvoluted spectrum to an SQLite database. The
    // MzmlMsSqlWriter consumes a raw MzmlMs (header + peaks) and reads the
    // activation's N/C ion types from the header, so an activation is attached
    // here. The text-peak-list input carries no activation information, so the
    // requested activation is used, falling back to HCD when it is unset
    // (the default "FILE" is not an activation name, so it is not looked up).
    if (para_ptr->isGeneSql()) {
      std::string activation_name = para_ptr->getActivation();
      if (activation_name == "FILE") {
        activation_name = "HCD";
      }
      ActivationPtr activation_ptr =
          ActivationBase::getActivationPtrByName(activation_name);
      if (activation_ptr == nullptr) {
        activation_ptr = ActivationBase::getActivationPtrByName("HCD");
      }
      header_ptr->setActivationPtr(activation_ptr);

      MzmlMsPtr raw_ms_ptr =
          std::make_shared<Ms<PeakPtr>>(header_ptr, peak_list);

      MzmlMsSqlWriterPtr sql_writer_ptr =
          std::make_shared<MzmlMsSqlWriter>(para_ptr->getSqlDb());
      sql_writer_ptr->writeMs2(raw_ms_ptr, result_envs);
      sql_writer_ptr->flush();
    }
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
}

}  // namespace

int process(const TopfdParaPtr& para_ptr,
            const std::vector<std::string>& spec_file_list) {
  // init data, envelope base, envcnn model, and ecscore model
  base_data::init(para_ptr->getResourceDir());
  EnvBase::initBase(para_ptr->getResourceDir());
  onnx_env_cnn::initModel(para_ptr->getResourceDir(), para_ptr->getThreadNum());

  for (const std::string& spec_file_name : spec_file_list) {
    std::cout << "Processing " << spec_file_name << " started." << std::endl;
    processOneFile(para_ptr, spec_file_name);
    std::cout << "Processing " << spec_file_name << " finished." << std::endl;
    std::cout << "Timestamp: " << time_util::getTimeStr() << std::endl;
  }

  base_data::release();
  std::cout << "TopFD single finished." << std::endl << std::flush;
  return 0;
}

}  // namespace topfd_single_process

}  // namespace toppic
