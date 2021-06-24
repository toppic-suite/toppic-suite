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

#ifndef TOPPIC_MERGE_FEATURE_PRSM_HPP_
#define TOPPIC_MERGE_FEATURE_PRSM_HPP_

#include "prsm/prsm_str.hpp"
#include "ms/feature/sample_feature.hpp"

namespace toppic {

class FeaturePrsm;
typedef std::shared_ptr<FeaturePrsm> FeaturePrsmPtr;

class FeaturePrsm : public SampleFeature {
 public:
  FeaturePrsm(std::string line);

  void addPrsmInfo(PrsmStrPtr prsm);

  std::string getProtName() {return prot_name_;}

  std::string getProtDesc() {return prot_desc_;}

  int getFirstResidue() {return first_residue_;}

  int getLastResidue() {return last_residue_;}

  std::string getProteoform() {return proteoform_;}

  int getMs2Id() {return ms2_id_;}

  double getPrecMass() {return prec_mass_;}

  double getAlignTimeBegin() {return align_time_begin_;}

  double getAlignTimeEnd() {return align_time_end_;}

  double getAlignTimeMiddle() {return (align_time_begin_ + align_time_end_)/2;}

  void setAlignTimeBegin(double time_begin) {align_time_begin_ = time_begin;}

  void setAlignTimeEnd(double time_end) {align_time_end_ = time_end;}

  static bool cmpMassInc(const FeaturePrsmPtr &a, const FeaturePrsmPtr &b) { 
    return a->getMonoMass() < b->getMonoMass();}

  static bool cmpInteDec(const FeaturePrsmPtr &a, const FeaturePrsmPtr &b) { 
    return a->getIntensity() > b->getIntensity();}

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
};

typedef std::vector<FeaturePrsmPtr> FeaturePrsmPtrVec;
typedef std::vector<FeaturePrsmPtrVec> FeaturePrsmPtrVec2D;

}

#endif
