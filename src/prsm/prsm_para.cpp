// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <boost/algorithm/string.hpp>

#include "base/activation_base.hpp"
#include "base/mass_constant.hpp"
#include "base/mod_util.hpp"
#include "base/prot_mod_base.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {

PrsmPara::PrsmPara(std::map<std::string, std::string> &arguments) {
  search_db_file_name_ = arguments["databaseFileName"];
  spec_file_name_ = arguments["spectrumFileName"];
  exe_dir_ = arguments["executiveDir"];
  errorTolerance_=std::stoi(arguments["errorTolerance"]);

  group_spec_num_ = std::stoi(arguments["groupSpectrumNumber"]);

  fix_mod_list_ = ModUtil::geneFixedModList(arguments["fixedMod"]);

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
  double IM = MassConstant::getIsotopeMass();
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

  sp_para_ptr_ = std::make_shared<SpPara>(min_peak_num, min_mass, extend_min_mass,
                                          ext_offsets, peak_tolerance_ptr, activation_ptr);
}

} /* namespace prot */
