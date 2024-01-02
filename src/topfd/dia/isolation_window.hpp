//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_TOPFD_DIA_ISOLATION_WINDOW_HPP_
#define TOPPIC_TOPFD_DIA_ISOLATION_WINDOW_HPP_

#include <memory>
#include <vector>
#include <set>

namespace toppic {

class IsolationWindow {
 public:
  IsolationWindow(double mz_bgn, double mz_end);

  void addSpecId(int id) {spec_id_set_.insert(id);}

  bool findSpecId(int id);

  bool isMatch(double mz_bgn, double mz_end);

 private:
  double mz_bgn_;
  double mz_end_;

  std::set<int> spec_id_set_;
};

typedef std::shared_ptr<IsolationWindow> IsolationWindowPtr;
typedef std::vector<IsolationWindowPtr> IsolationWindowPtrVec;

}

#endif
