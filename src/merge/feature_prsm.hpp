//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_MERGE_FEATURE_PRSM_HPP_
#define TOPPIC_MERGE_FEATURE_PRSM_HPP_

#include "prsm/prsm_str.hpp"

namespace toppic {

class FeaturePrsm;
typedef std::shared_ptr<FeaturePrsm> FeaturePrsmPtr;

class FeaturePrsm {
 public:
  FeaturePrsm(PrsmStrPtr prsm);

  int getSampleId() {return sample_id_;}

  int getProteoId() {return proteo_id_;}

  std::string getProtName() {return prot_name_;}

  std::string getProtDesc() {return prot_desc_;}

  int getFirstResidue() {return first_residue_;}

  int getLastResidue() {return last_residue_;}

  std::string getProteoform() {return proteoform_;}

  int getMs2Id() {return ms2_id_;}

  double getAlignMinTime() {return align_min_time_;}

  double getAlignMaxTime() {return align_max_time_;}

  double getAlignApexTime() {return align_apex_time_;}

  double getPrecMass() {return prec_mass_;}

  double getProteoInte() {return proteo_inte_;}

  double getMinTime() {return min_time_;}

  double getMaxTime() {return max_time_;}

  double getApexTime() {return apex_time_;}

  void setSampleId(int sample_id) {sample_id_ = sample_id;}

  void setAlignMinTime(double min_time) {align_min_time_ = min_time;}

  void setAlignMaxTime(double max_time) {align_max_time_ = max_time;}

  void setAlignApexTime(double apex_time) {align_apex_time_ = apex_time;}

  static bool cmpMassInc(const FeaturePrsmPtr &a, const FeaturePrsmPtr &b) { 
    return a->getPrecMass() < b->getPrecMass();}

  static bool cmpInteDec(const FeaturePrsmPtr &a, const FeaturePrsmPtr &b) { 
    return a->getProteoInte() > b->getProteoInte();}

  static bool cmpTimeInc(const FeaturePrsmPtr &a, const FeaturePrsmPtr &b) { 
    return a->getApexTime() < b->getApexTime();}

 private:
  std::string prot_name_;
  std::string prot_desc_;
  int first_residue_;
  int last_residue_;
  std::string proteoform_;
  int ms2_id_;
  double align_min_time_;
  double align_max_time_;
  double align_apex_time_;

  int sample_id_;

  // new variable needed from prsm_str
  double prec_mass_;
  int proteo_id_;
  double proteo_inte_;
  double max_time_;
  double min_time_;
  double apex_time_;
};

typedef std::vector<FeaturePrsmPtr> FeaturePrsmPtrVec;
typedef std::vector<FeaturePrsmPtrVec> FeaturePrsmPtrVec2D;

}

#endif
