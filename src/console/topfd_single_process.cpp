//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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
//
#include <iomanip>

#include "common/util/version.hpp"
#include "common/util/logger.hpp"
#include "common/util/time_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/base/base_data.hpp"
#include "common/base/mass_constant.hpp"
#include "ms/spec/msalign_frac_merge.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/feature/feature_merge.hpp"
#include "topfd/msreader/raw_ms_writer.hpp"
#include "topfd/deconv/deconv_process.hpp"
#include "topfd/feature_detect/feature_detect.hpp"

namespace toppic {

namespace topfd_single_process {

std::string geneArgumentStr(std::map<std::string, std::string> arguments, 
                            const std::string & prefix) {
  std::stringstream output;
  output << prefix << "TopFD " << version_number << std::endl;
  // TIME_STAMP_STR is replaced later
  output << prefix << "Timestamp: " << time_util::TIME_STAMP_STR << std::endl;
  output << prefix << "###################### Parameters ######################" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Data type: " << "centroid" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Precursor charge: " << arguments["precursorCharge"] << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Precursor monoisotopic mass: " << arguments["precursorMass"] << " Dalton" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Error tolerance: " << arguments["mzError"] << " m/z" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "MS/MS signal/noise ratio: " << arguments["msTwoSnRatio"] << std::endl;
  //output << prefix << std::setw(40) << std::left 
  //    << "Do final filtering: " << para_ptr->do_final_filtering_ << std::endl;
  output << prefix << "###################### Parameters ######################" << std::endl;
  return output.str();
}


PeakPtrVec readPeakFile(std::string file_name) {
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
    //std::cout << "mz " << mz << " intensity " << inte << std::endl;
    PeakPtr peak_ptr = std::make_shared<Peak>(mz, inte);
    peak_list.push_back(peak_ptr);
  }
  input.close();
  //std::cout << "peak list finished" << std::endl;
  return peak_list;
}

int processOneFile(std::map<std::string, std::string> arguments, 
                   const std::string &argument_str, 
                   const std::string &spec_file_name) {
  try {
    base_data::init(arguments["resourceDir"]);

    EnvParaPtr env_para_ptr = std::make_shared<EnvPara>();
    int charge = std::stoi(arguments["precursorCharge"]);
    env_para_ptr->max_charge_ = charge; 
    double mass = std::stod(arguments["precursorMass"]);
    env_para_ptr->max_mass_ = mass; 
    double tolerance = std::stod(arguments["mzError"]);
    env_para_ptr->setTolerance(tolerance);
    double sn = std::stod(arguments["msTwoSnRatio"]);
    env_para_ptr->ms_two_sn_ratio_ = sn;
    std::cout << "charge " << charge << " mass " << mass << " tolerance " << tolerance << " sn " << sn << std::endl;

    DpParaPtr dp_para_ptr = std::make_shared<DpPara>();
    DeconvOneSpPtr deconv_ptr = std::make_shared<DeconvOneSp>(env_para_ptr, dp_para_ptr);

    PeakPtrVec peak_list = readPeakFile(spec_file_name); 
    MsHeaderPtr header_ptr = std::make_shared<MsHeader>();
    header_ptr->setId(0);
    header_ptr->setFractionId(0);
    header_ptr->setScan(1);
    header_ptr->setPrecCharge(charge);
    double prec_mz = mass / charge + mass_constant::getProtonMass();
    header_ptr->setPrecMonoMz(prec_mz);
    header_ptr->setPrecSpMz(prec_mz);
    header_ptr->setMsLevel(2);
    header_ptr->setFileName(spec_file_name);
    RawMsPtr raw_ms_ptr = std::make_shared<Ms<PeakPtr>>(header_ptr, peak_list);
    std::cout <<"header finished" << std::endl;

    MatchEnvPtrVec result_envs; 
    if (peak_list.size() > 0) {
      deconv_ptr->setData(peak_list, mass, charge);
      deconv_ptr->run();
      result_envs = deconv_ptr->getResult();
    }
    DeconvMsPtr ms_ptr = match_env_util::getDeconvMsPtr(header_ptr, result_envs);

    std::cout <<"deconv finished" << std::endl;

    std::string ms2_msalign_name = file_util::basename(spec_file_name) + "_ms2.msalign";
    MsAlignWriterPtr ms2_writer_ptr = std::make_shared<MsAlignWriter>(ms2_msalign_name);
    ms2_writer_ptr->write(ms_ptr);
    ms2_writer_ptr->close();

    std::string spec_data_suffix = "_ms2_json";
    std::string json_dir =  file_util::basename(spec_file_name) + spec_data_suffix;
    file_util::createFolder(json_dir);

    std::string json_file_name = json_dir + file_util::getFileSeparator() + "spectrum" 
        + std::to_string(header_ptr->getId())
        + ".js";
    raw_ms_writer::write(json_file_name, raw_ms_ptr, result_envs);    

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopFD single finished." << std::endl << std::flush;
  return 0;
}

int process(std::map<std::string, std::string> arguments, 
            std::vector<std::string> spec_file_lst) {
  std::string argument_str = geneArgumentStr(arguments, "#");
  std::cout << argument_str;
  EnvBase::initBase(arguments["resourceDir"]);
  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    if (str_util::endsWith(spec_file_lst[k], "txt")) {
      int result = processOneFile(arguments, argument_str, spec_file_lst[k]);
      if (result != 0) {
        return 1;
      }
    }
  }
  return 0;
}


} // namespace topfd_process 

}  // namespace toppic

