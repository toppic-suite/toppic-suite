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


#ifndef PROT_SPEC_SP_PARA_HPP_
#define PROT_SPEC_SP_PARA_HPP_

#include <string>
#include <memory>
#include <vector>
#include <set>

#include "base/activation.hpp"
#include "xml/xml_dom_document.hpp"
#include "spec/peak_tolerance.hpp"

namespace toppic {

class SpPara {
 public:
  SpPara(int min_peak_num, double min_mass, double min_extend_mass,
         const std::vector<double> &ext_offsets,
         PeakTolerancePtr peak_tolerance_ptr,
         ActivationPtr activation_ptr,
         const std::set<std::string> & skip_list):
      min_peak_num_(min_peak_num),
      min_mass_(min_mass),
      extend_min_mass_(min_extend_mass),
      ext_offsets_(ext_offsets),
      peak_tolerance_ptr_(peak_tolerance_ptr),
      activation_ptr_(activation_ptr),
      skip_list_(skip_list) {}

  explicit SpPara(xercesc::DOMElement* element);

  double getMinMass() {return min_mass_;}

  double getExtendMinMass() {return extend_min_mass_;}

  const std::vector<double>& getExtendOffsets() {return ext_offsets_;}

  PeakTolerancePtr getPeakTolerancePtr() {return peak_tolerance_ptr_;}

  void setPeakTolerancePtr(PeakTolerancePtr peak_tolerance_ptr) {
    peak_tolerance_ptr_ = peak_tolerance_ptr;
  }

  ActivationPtr getActivationPtr() {return activation_ptr_;}

  void setActivationPtr(ActivationPtr activation_ptr) {activation_ptr_ = activation_ptr;}

  std::set<std::string> getSkipList() {return skip_list_;}

  int getMinPeakNum() {return min_peak_num_;}

  void setMinPeakNum(int min_peak_num) {min_peak_num_ = min_peak_num;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "sp_para";}

  // the 1 Da error in precursor mass used in zeroptm searching
  int prec_error_ = 1;

 private:
  int min_peak_num_;

  // if the mass if smaller than min_mass, the mass is removed.
  double min_mass_;
  // if the mass is smaller than extend_min_mass, the peak is not extended
  double extend_min_mass_;

  std::vector<double> ext_offsets_;

  PeakTolerancePtr peak_tolerance_ptr_;

  ActivationPtr activation_ptr_;

  std::set<std::string> skip_list_;
};

typedef std::shared_ptr<SpPara> SpParaPtr;

} /* namespace toppic */

#endif /* SP_PARA_HPP_ */
