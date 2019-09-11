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


#ifndef TOPPIC_DECONV_ENV_ENVELOPE_BASE_HPP_
#define TOPPIC_DECONV_ENV_ENVELOPE_BASE_HPP_

#include "deconv/env/envelope.hpp"

namespace toppic {

class EnvBase;

typedef std::shared_ptr<EnvBase> EnvBasePtr;

class EnvBase {
 public:
  EnvBase(std::string file_name_, int entry_num_, double mass_interval_);

  void initBaseMassIdx();

  std::vector<std::vector<double> > env_rescore_para_;

  static EnvelopePtr getStaticEnvByMonoMass(double mass);

  static EnvelopePtr getStaticEnvByBaseMass(double mass);

  static void initBase(const std::string &resouce_dir);

 private:
  // number of distribution entries 
  int entry_num_;
  // the mass interval between two neighboring entries 
  double mass_interval_;
  // the list of distribution envelopes 
  EnvelopePtrVec envs_;
  // mapping distribution envelopes to the mass value of base peak 
  std::vector<int> base_mass_idxes_;

  EnvelopePtr getEnvByMonoMass(double mass);

  EnvelopePtr getEnvByBaseMass(double mass);

  static EnvBasePtr env_base_ptr_;

  static double getDistrMassInterval() {return 10;}

  static int getDistrEntryNum() {return 11000;}

  static std::string getBaseDirName() {return "base_data";}

  static std::string getBaseFileName() {return "theo_patt.txt";}
};

}

#endif
