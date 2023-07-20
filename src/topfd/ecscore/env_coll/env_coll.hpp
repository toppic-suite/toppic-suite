//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_COLL_ENV_COLL_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_COLL_ENV_COLL_HPP

#include <memory>
#include <vector>

#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/envelope/seed_envelope.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"

namespace toppic {

class EnvColl {
 public:
  EnvColl(SeedEnvelopePtr send_ptr, EnvSetPtrVec &env_set_list, 
          int min_charge, int max_charge, 
          int start_spec_id, int end_spec_id);

  std::vector<double> compExpInteSumList(); 

  std::vector<int> getChargeList();

  void refineMonoMass();

  double getIntensity(double sn_ratio, double noise_inte);

  EnvSetPtr getSeedEnvSet();

  SeedEnvelopePtr getSeedPtr() { return seed_ptr_; }

  void setSeedPtr(SeedEnvelopePtr seed_ptr) { seed_ptr_ = seed_ptr; }

  EnvSetPtrVec getEnvSetList() { return env_set_list_; }

  void setEnvSetList(EnvSetPtrVec &env_set_list);  

  int getMinCharge() const { return min_charge_; }

  void setMinCharge(int min_charge) { min_charge_ = min_charge; }

  int getMaxCharge() const { return max_charge_; }

  void setMaxCharge(int max_charge) { max_charge_ = max_charge; }

  int getStartSpecId() const { return start_spec_id_; }

  void setStartSpecId(int start_spec_id) { start_spec_id_ = start_spec_id; }

  int getEndSpecId() const { return end_spec_id_; }

  void setEndSpecId(int end_spec_id) { end_spec_id_ = end_spec_id; }

  double getEcscore() const { return ecscore_; }

  void setEcscore(double ecscore) { ecscore_ = ecscore; }

  int getBaseSpecId() const { return seed_ptr_->getSpecId(); }

  double getMass() const { return seed_ptr_->getMonoNeutralMass(); }

  std::vector<double> getExpInteSumList() const { return exp_inte_sum_list_; }

  void setExpInteSumList(const std::vector<double> &exp_inte_sum_list) {
    exp_inte_sum_list_ = exp_inte_sum_list;}

  void removePeakData(MsMapPtr matrix_ptr);

 private:
  SeedEnvelopePtr seed_ptr_;
  EnvSetPtrVec env_set_list_;
  int min_charge_;
  int max_charge_;
  int start_spec_id_;
  int end_spec_id_;
  double ecscore_ = -1;
  std::vector<double> exp_inte_sum_list_;
};

typedef std::shared_ptr<EnvColl> EnvCollPtr;
typedef std::vector<EnvCollPtr> EnvCollPtrVec;


}  // namespace toppic

#endif //TOPPIC_ENV_COLLECTION_HPP
