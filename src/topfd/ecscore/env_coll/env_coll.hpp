//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include "common/xml/xml_dom_document.hpp"
#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/env/seed_env.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"

namespace toppic {

class EnvColl {
 public:
  EnvColl(SeedEnvPtr send_ptr, EnvSetPtrVec &env_set_list,
          int min_charge, int max_charge,
          int start_spec_id, int end_spec_id);

  SeedEnvPtr getSeedPtr() { return seed_ptr_; }

  EnvSetPtrVec getEnvSetList() { return env_set_list_; }

  int getMinCharge() const { return min_charge_; }

  int getMaxCharge() const { return max_charge_; }

  int getStartSpecId() const { return start_spec_id_; }

  int getEndSpecId() const { return end_spec_id_; }

  double getEcscore() const { return ecscore_; }

  void setEcscore(double ecscore) { ecscore_ = ecscore; }

  int getSeedSpecId() const { return seed_ptr_->getSpecId(); }

  double getMonoNeutralMass() const { return seed_ptr_->getMonoNeutralMass(); }

  double getIntensity();
  
  int countEnvNum();

  EnvSetPtr getSeedEnvSet();
  
  std::vector<int> getChargeList();

  void refineMonoMass();

  void removePeakData(MsMapPtr matrix_ptr);

  void mergeEnvSet(EnvSetPtr env_set_ptr);

  XmlDOMElement* toXmlElement(XmlDOMDocument* xml_doc);

 private:
  SeedEnvPtr seed_ptr_;
  EnvSetPtrVec env_set_list_;
  int min_charge_;
  int max_charge_;
  int start_spec_id_;
  int end_spec_id_;
  double ecscore_ = -1;
};

typedef std::shared_ptr<EnvColl> EnvCollPtr;
typedef std::vector<EnvCollPtr> EnvCollPtrVec;


}  // namespace toppic

#endif //TOPPIC_ENV_COLLECTION_HPP
