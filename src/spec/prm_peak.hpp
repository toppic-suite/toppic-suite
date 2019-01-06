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


#ifndef TOPPIC_SPEC_PRM_PEAK_HPP_
#define TOPPIC_SPEC_PRM_PEAK_HPP_

#include "spec/deconv_peak.hpp"
#include "spec/support_peak.hpp"
#include "spec/base_peak_type.hpp"
#include "spec/rm_break_type.hpp"

namespace toppic {

class PrmPeak;
typedef std::shared_ptr<PrmPeak> PrmPeakPtr;

class PrmPeak : public Peak {
 public:
  PrmPeak(int spec_id, DeconvPeakPtr base_peak_ptr,
          BasePeakTypePtr base_type,
          double mono_mass, double score,
          double strict_tolerance = 0.0,
          double n_strict_c_relax_tolerance = 0.0,
          double n_relax_c_strict_tolerance = 0.0);

  void addNghbEdge(DeconvPeakPtr deconv_peak_ptr, double offset,
                   SPTypePtr peak_type, double score);

  int getNeighborSize() {return neighbor_list_.size();}

  DeconvPeakPtr getBasePeakPtr() {return base_peak_ptr_;}

  double getMonoMass() {return mono_mass_;}

  void setMonoMass(double m);

  double getScore() {return score_;}

  double getStrictTolerance() {return strict_tolerance_;}

  BasePeakTypePtr getBaseTypePtr() {return base_type_;}

  double getNStrictCRelaxTolerance() {return n_strict_c_relax_tolerance_;}

  double getNRelaxCStrictTolerance() {return n_relax_c_strict_tolerance_;}

  int getSpectrumId() {return spec_id_;}

  int getPeakId() {return peak_id_;}

  RmBreakTypePtr getBreakType();

  void setStrictTolerance(double tolerance) {strict_tolerance_ = tolerance;}

  void setNStrictCRelacTolerance(double tolerance) {
    n_strict_c_relax_tolerance_ = tolerance;
  }

  void setNRelaxCStrictTolerance(double tolerance) {
    n_relax_c_strict_tolerance_ = tolerance;
  }

  void setPeakId(int peak_id) {peak_id_ = peak_id;}

  static bool cmpPosInc(const PrmPeakPtr &a, const PrmPeakPtr &b) {
    return a->getPosition() < b->getPosition();}

 private:
  int spec_id_;
  DeconvPeakPtr base_peak_ptr_;
  BasePeakTypePtr base_type_;
  int peak_id_;
  double mono_mass_;
  double score_;
  double strict_tolerance_;
  double n_strict_c_relax_tolerance_;
  double n_relax_c_strict_tolerance_;
  SupportPeakPtrVec neighbor_list_;
};

typedef std::vector<PrmPeakPtr> PrmPeakPtrVec;
typedef std::vector<PrmPeakPtrVec> PrmPeakPtrVec2D;


} /* namespace toppic */

#endif /* PRM_PEAK_HPP_ */
