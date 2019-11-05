//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef TOPPIC_PRSM_VIEW_ANNO_CLEAVAGE_HPP_
#define TOPPIC_PRSM_VIEW_ANNO_CLEAVAGE_HPP_

#include <string>
#include <vector>

#include "seq/proteoform.hpp"
#include "spec/extend_peak.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/prsm.hpp"

namespace toppic {

#define CLEAVAGE_TYPE_NORMAL "normal"
#define CLEAVAGE_TYPE_N_TRUNCATION "n_truncation"
#define CLEAVAGE_TYPE_C_TRUNCATION "c_truncation"
#define CLEAVAGE_TYPE_SEQ_START "seq_start"
#define CLEAVAGE_TYPE_SEQ_END "seq_end"

class AnnoCleavage {
 public:
  AnnoCleavage(int pos, const PeakIonPairPtrVec &pairs, bool exist_n_ion, bool exist_c_ion):
      pos_(pos),
      pairs_(pairs),
      exist_n_ion_(exist_n_ion),
      exist_c_ion_(exist_c_ion),
      is_unexpected_change_(false),
      unexpected_change_color_(0),
      type_(CLEAVAGE_TYPE_NORMAL) {}

  void setPairs(PeakIonPairPtrVec pairs) {pairs_ = pairs;}

  void setExistNIon(bool n) {exist_n_ion_ = n;}

  void setExistCIon(bool c) {exist_c_ion_ = c;}

  void setType(const std::string &type) {type_ = type;}

  void setUnexpectedChange(bool u) {is_unexpected_change_ = u;}

  void setUnexpectedChangeColor(int color) {unexpected_change_color_ = color;}

  std::string getType() {return type_;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
  int pos_;

  PeakIonPairPtrVec pairs_;

  bool exist_n_ion_;

  bool exist_c_ion_;

  bool is_unexpected_change_;

  int unexpected_change_color_;

  std::string type_;
};

typedef std::shared_ptr<AnnoCleavage> AnnoCleavagePtr;
typedef std::vector<AnnoCleavagePtr> AnnoCleavagePtrVec;

AnnoCleavagePtrVec getProteoCleavage(PrsmPtr prsm_ptr, double min_mass);
} /* namespace toppic */

#endif /* TOPPIC_PRSM_VIEW_ANNO_CLEAVAGE_HPP_ */
