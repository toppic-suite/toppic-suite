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


#include <boost/algorithm/string.hpp>

#include "base/activation_base.hpp"
#include "base/mass_constant.hpp"
#include "base/mod_util.hpp"
#include "base/prot_mod_base.hpp"
#include "prsm/prsm_para.hpp"

namespace toppic {

PrsmPara::PrsmPara(std::map<std::string, std::string> &arguments) {
  search_db_file_name_ = arguments["databaseFileName"];

  spec_file_name_ = arguments["spectrumFileName"];

  resource_dir_ = arguments["resourceDir"];

  errorTolerance_=std::stoi(arguments["errorTolerance"]);

  group_spec_num_ = std::stoi(arguments["groupSpectrumNumber"]);

  fix_mod_list_ = mod_util::geneFixedModList(arguments["fixedMod"]);

  std::string prot_mod_str = arguments["allowProtMod"];
  std::vector<std::string> strs;
  boost::split(strs, prot_mod_str, boost::is_any_of(","));
  for (size_t i = 0; i < strs.size(); i++) {
    ProtModPtrVec mods = ProtModBase::getProtModPtrByType(strs[i]);
    LOG_DEBUG("prot mod type " << strs[i] << " num " << mods.size());
    prot_mod_list_.insert(prot_mod_list_.end(), mods.begin(), mods.end());
  }

  std::string activation_name = arguments["activation"];
  ActivationPtr activation_ptr 
      = ActivationBase::getActivationPtrByName(activation_name);

  double ppo = std::stod(arguments["errorTolerance"])*0.000001;
  bool use_min_tolerance = true;
  double min_tolerance = 0.01;
  PeakTolerancePtr peak_tolerance_ptr
      = std::make_shared<PeakTolerance>(ppo, use_min_tolerance, min_tolerance);

  // extend sp parameter 
  double IM = mass_constant::getIsotopeMass();
  // the set of offsets used to expand the monoisotopic mass list 
  std::vector<double> ext_offsets {{0, -IM, IM}};
  double extend_min_mass = 5000;

  int min_peak_num = 10;
  double min_mass = 50.0;

  if (arguments["residueModFileName"] != "") {
    localization_ = true;
  } else {
    localization_ = false;
  }

  std::set<std::string> skip_list;

  if (arguments["skipList"] != "") {
    std::ifstream infile(arguments["skipList"]);
    std::string line;
    while (std::getline(infile, line)) {
      std::vector<std::string> strs; 
      boost::split(strs, line, boost::is_any_of(" "));
      for (size_t i = 0; i < strs.size(); i++) {
        skip_list.insert(strs[i]);
      }
    }
    infile.close();
  }

  sp_para_ptr_ = std::make_shared<SpPara>(min_peak_num, min_mass, extend_min_mass,
                                          ext_offsets, peak_tolerance_ptr,
                                          activation_ptr, skip_list);
}

} /* namespace toppic */
