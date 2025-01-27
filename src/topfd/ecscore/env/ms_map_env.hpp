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

#ifndef TOPPIC_ECSCORE_ENV_MS_MAP_ENV_HPP
#define TOPPIC_ECSCORE_ENV_MS_MAP_ENV_HPP

#include <vector>

#include "common/xml/xml_dom_document.hpp"
#include "ms/msmap/ms_map_peak.hpp"
#include "topfd/ecscore/env/seed_env.hpp"

namespace toppic {

class MsMapEnv {
 public:
  MsMapEnv(int spec_id, MsMapPeakPtrVec peak_list);

  std::vector<double> getInteList();

  std::vector<double> getMzList();

  int getPeakNum() { return peak_list_.size(); }

  double getInteSum(); 

  MsMapPeakPtr getPeakPtr(int idx) { return peak_list_[idx]; }

  int getSpecId() const { return spec_id_; }

  void setSpecId(int spec_id) { spec_id_ = spec_id; }

  MsMapPeakPtrVec getMsMapPeakList() { return peak_list_; }

  void setPeakPtr(int idx, MsMapPeakPtr peak_ptr) { peak_list_[idx] = peak_ptr; }

  int getTopThreeMatchNum(int ref_idx);

  double compTopThreeInteSum(int ref_idx);

  void removeLowIntePeaks(SeedEnvPtr seed_ptr, double ratio, double min_inte);

  void appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

 private:
  int spec_id_;
  MsMapPeakPtrVec peak_list_;
};

typedef std::shared_ptr<MsMapEnv> MsMapEnvPtr;
typedef std::vector<MsMapEnvPtr> MsMapEnvPtrVec;

}

#endif
