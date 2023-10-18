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

#ifndef TOPPIC_TOPFD_ENV_ENVELOPE_BASE_HPP_
#define TOPPIC_TOPFD_ENV_ENVELOPE_BASE_HPP_

#include "ms/env/env.hpp"

namespace toppic {

class EnvBase;

typedef std::shared_ptr<EnvBase> EnvBasePtr;

// A list of envelopes read from a file theo_patt.txt
class EnvBase {
 public:
  EnvBase(std::string file_name_, int entry_num_, double mass_interval_);

  static void initBase(const std::string &resouce_dir);

  // All public functions are static 
  static EnvPtr getEnvByMonoMass(double mass);

  static EnvPtr getEnvByMonoMass(double mass, int charge);

  static EnvPtr getEnvByRefMass(double mass);

  static double convertMonoMassToAvgMass(double mass);

  static double convertMonoMassToRefMass(double mass);

  static double convertRefMassToMonoMass(double mass); 

 private:
  // number of distribution entries 
  int entry_num_;
  // the mass interval between two neighboring entries 
  double mass_interval_;
  // the list of distribution envelopes 
  EnvPtrVec envs_;
  // mapping distribution envelopes to the mass value of reference peak (highest intensity peak) 
  std::vector<int> ref_mass_idxes_;

  static EnvBasePtr env_base_ptr_;

  static double getDistrMassInterval() {return 10;}

  static int getDistrEntryNum() {return 11000;}

  static std::string getBaseDirName() {return "base_data";}

  static std::string getBaseFileName() {return "theo_patt.txt";}

  void initRefMassIdx();
  
  // this is a private class method
  EnvPtr getBaseEnvByMonoMass(double mass);

  // get a envelope using the mass of the highest peak
  EnvPtr getBaseEnvByRefMass(double mass);
};

}

#endif
