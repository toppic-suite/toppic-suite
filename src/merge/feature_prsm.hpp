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

  double getPrecMass() {return prec_mass_;}

  double getAlignTimeBegin() {return align_time_begin_;}

  double getAlignTimeEnd() {return align_time_end_;}

  //double getAlignTimeMiddle() {return (align_time_begin_ + align_time_end_)/2;}
  double getAlignTimeApex() {return align_time_apex_;}

  double getMonoMass() {return mono_mass_;}

  double getIntensity() {return inte_;}

  double getTimeBegin() {return time_begin_;}

  double getTimeEnd() {return time_end_;}

  double getTimeApex() {return time_apex_;}

  void setSampleId(int sample_id) {sample_id_ = sample_id;}

  void setAlignTimeBegin(double time_begin) {align_time_begin_ = time_begin;}

  void setAlignTimeEnd(double time_end) {align_time_end_ = time_end;}

  void setAlignTimeApex(double time_apex) {align_time_apex_ = time_apex;}

  static bool cmpMassInc(const FeaturePrsmPtr &a, const FeaturePrsmPtr &b) { 
    return a->getMonoMass() < b->getMonoMass();}

  static bool cmpInteDec(const FeaturePrsmPtr &a, const FeaturePrsmPtr &b) { 
    return a->getIntensity() > b->getIntensity();}

  static bool cmpTimeInc(const FeaturePrsmPtr &a, const FeaturePrsmPtr &b) { 
    return a->getTimeApex() < b->getTimeApex();}

 private:
  std::string prot_name_;
  std::string prot_desc_;
  int first_residue_;
  int last_residue_;
  std::string proteoform_;
  int ms2_id_;
  double prec_mass_;
  double align_time_begin_;
  double align_time_end_;
  double align_time_apex_;

  int sample_id_;

  // new variable needed from prsm_str
  double mono_mass_;
  double inte_;
  double time_end_;
  double time_begin_;
  double time_apex_;
  int proteo_id_;
};

typedef std::vector<FeaturePrsmPtr> FeaturePrsmPtrVec;
typedef std::vector<FeaturePrsmPtrVec> FeaturePrsmPtrVec2D;

}

#endif
