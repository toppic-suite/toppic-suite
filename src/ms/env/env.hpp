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

#ifndef TOPPIC_TOPFD_ENV_ENV_HPP_
#define TOPPIC_TOPFD_ENV_ENV_HPP_

#include "ms/spec/peak_util.hpp"
#include "ms/spec/env_peak.hpp"

namespace toppic {

class Env;

typedef std::shared_ptr<Env> EnvPtr;

class Env {
 public:
  Env() {}

  Env(Env &env);

  //used by env_base to create envelopes
  Env(int num, std::vector<std::string> &line_List);

  Env(int refer_idx, int charge, double mono_mz, EnvPeakPtrVec &peaks);

  EnvPtr convertToTheo(double mass_diff, int new_charge);

  EnvPtr distrToTheoRef(double new_ref_mz, int new_charge);

  EnvPtr distrToTheoMono(double new_mono_mz, int new_charge);

  void changeIntensity(double ratio);

  void changeToAbsInte(double absolute_intensity);

  void changeMz(double shift);

  EnvPtr getSubEnv(int n_back, int n_forw);

  EnvPtr getSubEnv(double min_inte);

  EnvPtr addZero(int num);

  EnvPtr getSubEnv(double percent_bound, double absolute_min_inte,
                   int max_back_peak_num, int max_forw_peak_num);

  std::vector<int> calcBound(double percent_bound, double absolute_min_inte,
                             int max_back_peak_num, int max_forw_peak_num);

  void shift(int shift);

  double getMinMz() {return peak_ptr_list_[0]->getPosition();}

  double getMaxMz() {return peak_ptr_list_[peak_ptr_list_.size() - 1]->getPosition();}

  double getMonoMz() {return mono_mz_;}

  double getMonoNeutralMass() {return peak_util::compPeakNeutralMass(mono_mz_, charge_);}

  // get the m/z difference between mono_mz and reference peak 
  double getMonoReferDistance() {return peak_ptr_list_[refer_idx_]->getPosition() - mono_mz_;}

  double getAvgMz();

  double getAvgNeutralMass() {return peak_util::compPeakNeutralMass(getAvgMz(), charge_);}

  double getReferMz() {return peak_ptr_list_[refer_idx_]->getPosition();}

  double getReferNeutralMass() {return peak_util::compPeakNeutralMass(getReferMz(), charge_);}

  double getMz(int i) {return peak_ptr_list_[i]->getPosition();}

  std::vector<double> getMzList();

  int getReferIdx() {return refer_idx_;}

  double getInte(int i) {return peak_ptr_list_[i]->getIntensity();}

  double getReferInte() {return peak_ptr_list_[refer_idx_]->getIntensity();}

  std::vector<double> getInteList();

  double compInteSum();

  int getCharge() {return charge_;}

  int getPeakNum() {return peak_ptr_list_.size();}

  EnvPeakPtr getPeakPtr(int i) {return peak_ptr_list_[i];}

  EnvPeakPtrVec getPeakPtrList() {return peak_ptr_list_;}

  void setPeakPtrList(EnvPeakPtrVec peak_ptr_list) { peak_ptr_list_ = peak_ptr_list;}

  void setInte(int i, double intensity) {peak_ptr_list_[i]->setIntensity(intensity);}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "envelope";}

 protected:
  // the peak index with the highest intensity
  int refer_idx_;
  // Charge of the envolope 
  int charge_;
  // Theoretical m/z value of monoisotopic ion 
  double mono_mz_;
  // peak list
  EnvPeakPtrVec peak_ptr_list_;

  //used in finding the reference index
  int getHighestPeakIdx();
};

typedef std::vector<EnvPtr> EnvPtrVec;

}

#endif
