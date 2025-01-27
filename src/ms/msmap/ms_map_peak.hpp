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

#ifndef TOPPIC_MS_MS_MAP_MS_MAP_PEAK_HPP
#define TOPPIC_MS_MS_MAP_MS_MAP_PEAK_HPP

#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "ms/spec/peak.hpp"

namespace toppic {

class MsMapPeak : public Peak {
 public:
  MsMapPeak(PeakPtr peak);

  std::string getString();

  double getOriInte() const { return ori_inte_; }

  void setOriInte(double ori_inte) { ori_inte_ = ori_inte; }

  bool getNeighbor() const { return neighbor_; }

  void setNeighbor(bool neighbor) { neighbor_ = neighbor; }

  void appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

 private:
  double ori_inte_;
  bool neighbor_;
};

typedef std::shared_ptr<MsMapPeak> MsMapPeakPtr;
typedef std::vector<MsMapPeakPtr> MsMapPeakPtrVec;
typedef std::vector<MsMapPeakPtrVec> MsMapPeakPtr2D;

}
#endif
