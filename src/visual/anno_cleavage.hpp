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


#ifndef TOPPIC_VISUAL_ANNO_CLEAVAGE_HPP_
#define TOPPIC_VISUAL_ANNO_CLEAVAGE_HPP_

#include <string>
#include <vector>

#include "prsm/peak_ion_pair.hpp"
#include "prsm/prsm.hpp"

namespace toppic {

class AnnoCleavage;
typedef std::shared_ptr<AnnoCleavage> AnnoCleavagePtr;
typedef std::vector<AnnoCleavagePtr> AnnoCleavagePtrVec;

class AnnoCleavage {
 public:
  AnnoCleavage(int pos, const PeakIonPairPtrVec &pairs, 
               bool exist_n_ion, bool exist_c_ion);

  void setPairs(PeakIonPairPtrVec pairs) {pairs_ = pairs;}

  void setExistNIon(bool n) {exist_n_ion_ = n;}

  void setExistCIon(bool c) {exist_c_ion_ = c;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static AnnoCleavagePtrVec getProteoCleavage(PrsmPtr prsm_ptr, double min_mass);

 private:
  int pos_;

  PeakIonPairPtrVec pairs_;

  bool exist_n_ion_;

  bool exist_c_ion_;

};


} /* namespace toppic */

#endif /* TOPPIC_VISUAL_ANNO_CLEAVAGE_HPP_ */
