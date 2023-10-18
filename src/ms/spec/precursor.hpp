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

#ifndef TOPPIC_MS_SPEC_PRECURSOR_HPP_
#define TOPPIC_MS_SPEC_PRECURSOR_HPP_

#include <vector>
#include <memory>

#include "common/xml/xml_dom_document.hpp"

namespace toppic {

class Precursor;
typedef std::shared_ptr<Precursor> PrecursorPtr;

class Precursor {
 public:
  Precursor(int id, int feat_id, 
            double mono_mz, int charge, double inte);

  Precursor(XmlDOMElement* element);

  int getPrecId() {return prec_id_;}

  int getFeatureId() {return feat_id_;}

  double getMonoMz();

  double getAdjustedMonoMz() {return adjusted_mono_mz_;}

  int getCharge() {return charge_;}

  double getInte() {return inte_;}

  double getMonoMass();

  double getMonoMassMinusWater();

  double getErrorTolerance(double ppo) {return getMonoMass() * ppo;}

  std::pair<int,int> getMonoMassMinusWaterError(double ppo, double scale);

  void setPrecId(int prec_id) {prec_id_ = prec_id;}

  void setMonoMz(double mono_mz) {mono_mz_ = mono_mz;}

  void setAdjustedMonoMz(double adjusted_mono_mz) {adjusted_mono_mz_ = adjusted_mono_mz;}

  void setCharge(int charge) {charge_ = charge;}

  void setInte(double inte) {inte_ = inte;}

  XmlDOMElement* getPrecursorXml(XmlDOMDocument* xml_doc);

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static std::string getXmlElementName() {return "precursor";}

  static bool cmpInteDec(const PrecursorPtr &a, const PrecursorPtr &b);

 private:
  int prec_id_ = -1;
  // fraction feature id
  int feat_id_ = -1;
  // computed monoisotopic precursor m/z value 
  double mono_mz_ = -1;
  // adjusted mono_mz after proteoform identification
  double adjusted_mono_mz_ = -1;
  // precursor charge state  
  int charge_ = -1;
  // precursor intensity 
  double inte_ = 0;
};

typedef std::vector<PrecursorPtr> PrecursorPtrVec;

}

#endif
