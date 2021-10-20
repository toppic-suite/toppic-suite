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

#ifndef TOPPIC_SEQ_LOCAL_ANNO_HPP_
#define TOPPIC_SEQ_LOCAL_ANNO_HPP_

#include "common/base/ptm.hpp"
#include "common/xml/xml_dom_document.hpp"

namespace toppic {

// Contain PTM localization results of the MIScore method 
class LocalAnno {
 public:
  explicit LocalAnno(XmlDOMElement* element);

  LocalAnno(int left_pos, int right_pos, double conf,
            const std::vector<double> & scr_vec,
            double raw_scr, PtmPtr p);

  int getLeftBpPos() {return left_pos_;}

  int getRightBpPos() {return right_pos_;}
  // get the confidence score
  double getConf() {return conf_;}
  // get the score vector containing score for each site
  std::vector<double> getScrVec() {return scr_vec_;}

  double getScr();

  PtmPtr getPtmPtr() {return ptm_ptr_;}

  // raw score is not used in toppic and topmg
  double getRawScr() {return raw_scr_;}

  void setRawScr(double s) {raw_scr_ = s;}

  double getMassShift() {return mass_shift_;}

  void setMassShift(double m) {mass_shift_ = m;}

  static std::string getXmlElementName() {return "localization_annotation";}

  void appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

 private:
  // left and right position filtered by threshold
  int left_pos_, right_pos_;
  // confidence score
  double conf_;
  // score for each possible site
  std::vector<double> scr_vec_;
  // raw score
  double raw_scr_;
  // PTM
  PtmPtr ptm_ptr_;
  // mass shift
  double mass_shift_;
};

typedef std::shared_ptr<LocalAnno> LocalAnnoPtr;
typedef std::vector<LocalAnnoPtr> LocalAnnoPtrVec;

}  // namespace toppic

#endif
