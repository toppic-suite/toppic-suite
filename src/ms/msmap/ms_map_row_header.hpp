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

#ifndef TOPPIC_MS_MSMAP_MS_MAP_ROW_HEADER_HPP_
#define TOPPIC_MS_MSMAP_MS_MAP_ROW_HEADER_HPP_

#include <memory>
#include <vector>

namespace toppic {

class MsMapRowHeader {
 public:
  MsMapRowHeader(int spec_id, int scan_num, double rt);

  int getSpecId() const { return spec_id_; }
  void setSpecId(int specId) { spec_id_ = specId; }

  int getScanNum() const { return scan_num_; }
  void setScanNum(int scanNum) { scan_num_ = scanNum; }

  double getRt() const { return rt_; }
  void setRt(double rt) { rt_ = rt; }

  double getBaseInte() const { return base_inte_; }
  void setBaseInte(double base_inte) { base_inte_ = base_inte; }

 private:
  int spec_id_;
  int scan_num_;
  double rt_;
  double base_inte_;
};

using MsMapRowHeaderPtr = std::shared_ptr<MsMapRowHeader>;
using MsMapRowHeaderPtrVec = std::vector<MsMapRowHeaderPtr>;

}  // namespace toppic

#endif
