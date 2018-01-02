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


#ifndef PROT_FEATURE_ENVELOPE_HPP_
#define PROT_FEATURE_ENVELOPE_HPP_

#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include "spec/env_peak.hpp"

namespace prot {

class Envelope;

typedef std::shared_ptr<Envelope> EnvelopePtr;

class Envelope {
 public:
  Envelope() {}

  Envelope(Envelope &env); 

  Envelope(int num, std::vector<std::string> &line_List);

  Envelope(int refer_idx, int charge, double mono_mz, EnvPeakPtrVec &peaks);

  EnvelopePtr convertToTheo(double mass_diff, int new_charge);

  EnvelopePtr distrToTheoBase(double new_base_mz, int new_charge);

  EnvelopePtr distrToTheoMono(double new_mono_mz, int new_charge);

  void changeIntensity(double ratio);

  void changeToAbsInte(double absolute_intensity);

  void changeMz(double shift);

  EnvelopePtr getSubEnv(int n_back, int n_forw);

  EnvelopePtr addZero(int num);

  EnvelopePtr getSubEnv(double percent_bound, double absolute_min_inte,
                        int max_back_peak_num, int max_forw_peak_num);

  std::vector<int> calcBound(double percent_bound, double absolute_min_inte,
                             int max_back_peak_num, int max_forw_peak_num);

  void shift(int shift);

  double compIntensitySum();

  double getAvgMz();

  double getAvgMass();

  int getHighestPeakIdx();

  std::vector<double> getIntensities();

  int getCharge() {return charge_;}

  int getLabel(int i) {return (int)std::round((peaks_[i]->getPosition() - mono_mz_) * charge_);}

  double getIntensity(int i) {return peaks_[i]->getIntensity();}

  double getMonoMass() {return Peak::compPeakMass(mono_mz_, charge_);}

  double getMonoMz() {return mono_mz_;}

  // get the m/z difference between mono_mz and reference peak 
  double getMonoReferDistance() {return peaks_[refer_idx_]->getPosition() - mono_mz_;}

  double getMz(int i) {return peaks_[i]->getPosition();}

  int getPeakNum() {return peaks_.size();}

  EnvPeakPtr getPeakPtr(int i) {return peaks_[i];}

  int getReferIdx() {return refer_idx_;}

  double getReferIntensity() {return peaks_[refer_idx_]->getIntensity();}

  double getReferMz() {return peaks_[refer_idx_]->getPosition();}

  void setIntensity(int i, double intensity) {peaks_[i]->setIntensity(intensity);}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "envelope";}

 protected:

  int refer_idx_;
  // Charge of the envolope 
  int charge_;
  // Theoretical m/z value of monoisotopic ion 
  double mono_mz_;
  // peak list
  EnvPeakPtrVec peaks_;
};

typedef std::vector<EnvelopePtr> EnvelopePtrVec;

}

#endif
