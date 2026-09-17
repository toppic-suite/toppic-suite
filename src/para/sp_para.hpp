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

#ifndef TOPPIC_PARA_SP_PARA_HPP_
#define TOPPIC_PARA_SP_PARA_HPP_

#include <memory>
#include <string>
#include <vector>

#include "common/base/activation.hpp"
#include "common/xml/xml_dom_element.hpp"
#include "para/peak_tolerance.hpp"

namespace toppic {

class XmlDOMDocument;

class SpPara {
 public:
  SpPara(const std::string& activation_name, double n_term_mod_mass,
         double ppm);

  explicit SpPara(XmlDOMElement element);

  double getMinMass() const { return min_mass_; }

  double getExtendMinMass() const { return extend_min_mass_; }

  const std::vector<double>& getFragExtOffsets() const {
    return frag_ext_offsets_;
  }

  const std::vector<double>& getZeroShiftSearchPrecErrorVec() const {
    return zero_shift_search_prec_error_vec_;
  }

  const std::vector<double>& getVarPtmSearchPrecErrorVec() const {
    return var_ptm_search_prec_error_vec_;
  }

  const std::vector<double>& getOneShiftSearchPrecErrorVec() const {
    return one_shift_search_prec_error_vec_;
  }

  const std::vector<double>& getMultiShiftSearchPrecErrorVec() const {
    return multi_shift_search_prec_error_vec_;
  }

  PeakTolerancePtr getPeakTolerancePtr() const { return peak_tolerance_ptr_; }

  ActivationPtr getActivationPtr() const { return activation_ptr_; }

  int getMinPeakNum() const { return min_peak_num_; }

  double getNTermLabelMass() const { return n_term_label_mass_; }

  double getEnvCnnCutoff() const { return env_cnn_cutoff_; }

  void setEnvCnnCutoff(double cutoff) { env_cnn_cutoff_ = cutoff; }

  // Allow (true) or disallow (false) a +-1 Da precursor mass error in the
  // zero-shift, one-shift and variable-PTM filtering and search stages.
  void setZeroShiftPrecOneDaError(bool allow);
  void setOneShiftPrecOneDaError(bool allow);
  void setVarPtmPrecOneDaError(bool allow);

  // Extend (true) each fragment mass above extend_min_mass_ by a +-1 Da
  // isotope error, i.e. frag_ext_offsets_ = {0, -1 Da, +1 Da}; false keeps
  // only the mass itself ({0}).
  void setFragMassExtension(bool extend);

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement parent) const;

  static std::string getXmlElementName() { return "sp_para"; }

  static int getMaxSpecNumPerFile() { return 10000000; }

  static int getMaxFeatureNumPerFile() { return 10000000; }

 private:
  int min_peak_num_ = 10;
  // if the mass if smaller than min_mass, the mass is removed.
  double min_mass_ = 50.0;
  // if the mass is smaller than extend_min_mass, the peak is not extended
  double extend_min_mass_ = 5000;

  std::vector<double> frag_ext_offsets_;

  // n_term_label_mass is for iTRAQ or TMT labeling
  double n_term_label_mass_ = 0.0;

  // envelopes with an EnvCNN score below env_cnn_cutoff are removed
  double env_cnn_cutoff_ = 0.2;

  ActivationPtr activation_ptr_;

  // the precursor mass offsets tried by each filtering/search stage: {0}
  // (exact) or {0, -1 Da, +1 Da} (a +-1 Da precursor mass error is allowed)
  std::vector<double> zero_shift_search_prec_error_vec_;

  std::vector<double> one_shift_search_prec_error_vec_;

  std::vector<double> multi_shift_search_prec_error_vec_;

  std::vector<double> var_ptm_search_prec_error_vec_;

  PeakTolerancePtr peak_tolerance_ptr_;
};

using SpParaPtr = std::shared_ptr<SpPara>;
using SpParaPtrVec = std::vector<SpParaPtr>;

}  // namespace toppic

#endif
