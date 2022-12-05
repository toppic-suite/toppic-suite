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

#ifndef TOPPIC_PARA_SP_PARA_HPP_
#define TOPPIC_PARA_SP_PARA_HPP_

#include <string>
#include <memory>
#include <vector>

#include "common/xml/xml_dom_element.hpp"
#include "common/base/activation.hpp"
#include "para/peak_tolerance.hpp"

namespace toppic {

class XmlDOMDocument;

class SpPara {
 public:

  SpPara(std::string activation_name, double n_term_mod_mass, 
         double ppm);

  explicit SpPara(xercesc::DOMElement* element);

  double getMinMass() {return min_mass_;}

  double getExtendMinMass() {return extend_min_mass_;}

  const std::vector<double>& getExtendOffsets() {return ext_offsets_;}

  const std::vector<double>& getZeroShiftSearchPrecErrorVec() {
    return zero_shift_search_prec_error_vec_;}

  const std::vector<double>& getVarPtmSearchPrecErrorVec() {
    return var_ptm_search_prec_error_vec_;}

  const std::vector<double>& getOneShiftSearchPrecErrorVec() {
    return one_shift_search_prec_error_vec_;}

  const std::vector<double>& getMultiShiftSearchPrecErrorVec() {
    return multi_shift_search_prec_error_vec_;}

  PeakTolerancePtr getPeakTolerancePtr() {return peak_tolerance_ptr_;}

  ActivationPtr getActivationPtr() {return activation_ptr_;}

  int getMinPeakNum() {return min_peak_num_;}

  double getNTermLabelMass() {return n_term_label_mass_;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "sp_para";}

  static int getMaxSpecNumPerFile() {return 1000000;}

  static int getMaxFeatureNumPerFile() {return 1000000;}

 private:
  int min_peak_num_ = 10;
  // if the mass if smaller than min_mass, the mass is removed.
  double min_mass_ = 50.0;
  // if the mass is smaller than extend_min_mass, the peak is not extended
  double extend_min_mass_ = 5000;

  std::vector<double> ext_offsets_;

  // n_term_label_mass is for iTRAQ or TMT labeling
  double n_term_label_mass_ = 0.0;

  ActivationPtr activation_ptr_;

  // the 1 Da error in precursor mass used in zeroptm filtering
  std::vector<double> zero_shift_search_prec_error_vec_;

  std::vector<double> one_shift_search_prec_error_vec_;

  std::vector<double> multi_shift_search_prec_error_vec_;

  std::vector<double> var_ptm_search_prec_error_vec_;

  PeakTolerancePtr peak_tolerance_ptr_;
};

typedef std::shared_ptr<SpPara> SpParaPtr;
typedef std::vector<SpParaPtr> SpParaPtrVec;

} /* namespace toppic */

#endif /* SP_PARA_HPP_ */
