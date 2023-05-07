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

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_SET_XIC_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_SET_XIC_HPP

#include <memory>
#include <vector>

namespace toppic {

class Xic {
 public:
  Xic();

  Xic(int start_spec_id, int base_spec_id, 
      std::vector<double> &inte_list, 
      std::vector<double> &env_inte_list);

  Xic(int start_spec_id, int base_spec_id, 
      std::vector<double> &inte_list, 
      std::vector<double> &smooth_inte_list,
      std::vector<double> &env_inte_list);

  void moving_avg(int n);

  int getStartSpecId() { return start_spec_id_; }

  void setStartSpecId(int startSpecId) { start_spec_id_ = startSpecId; }

  int getBaseSpecId() { return base_spec_id_; }

  void setBaseSpecId(int baseSpecId) { base_spec_id_ = baseSpecId; }

  std::vector<double> getInteList() { return inte_list_; }

  void setInteList(std::vector<double> inteList); 

  double get_inte_list_sum() {return inte_sum_; }

  std::vector<double> getSmoothedInteList() { return smoothed_inte_list_; }

  void setSmoothedInteList(std::vector<double> smoothedInteList) { smoothed_inte_list_ = smoothedInteList; }

  std::vector<double> getEnvInteList() { return env_inte_list_; }

  void setEnvInteList(std::vector<double> envInteList) { env_inte_list_ = envInteList; }

 private:
  int start_spec_id_;
  int base_spec_id_;
  std::vector<double> inte_list_;
  std::vector<double> smoothed_inte_list_;
  std::vector<double> env_inte_list_;
  double inte_sum_;
};

typedef std::shared_ptr<Xic> XicPtr;
typedef std::vector<XicPtr> XicPtrVec;

}

#endif //TOPPIC_XIC_HPP
