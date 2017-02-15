// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_SPEC_PRM_PEAK_HPP_
#define PROT_SPEC_PRM_PEAK_HPP_

#include <memory>
#include <vector>

#include "spec/deconv_peak.hpp"
#include "spec/support_peak.hpp"
#include "spec/base_peak_type.hpp"
#include "spec/rm_break_type.hpp"

namespace prot {

class PrmPeak;
typedef std::shared_ptr<PrmPeak> PrmPeakPtr;

class PrmPeak : public Peak {
 public:
  PrmPeak(int spec_id, DeconvPeakPtr base_peak_ptr,
          BasePeakTypePtr base_type, 
          double mono_mass, double score);

  void addNghbEdge(DeconvPeakPtr deconv_peak_ptr, double offset,
                   SPTypePtr peak_type, double score);

  int getNeighborSize(){return neighbor_list_.size();}

  DeconvPeakPtr getBasePeakPtr(){return base_peak_ptr_;}

  double getMonoMass(){return mono_mass_;}

  double getScore(){return score_;}

  double getStrictTolerance(){return strict_tolerance_;}

  BasePeakTypePtr getBaseTypePtr(){return base_type_;}

  double getNStrictCRelaxTolerance(){return n_strict_c_relax_tolerance_;}

  double getNRelaxCStrictTolerance(){return n_relax_c_strict_tolerance_;}

  int getSpecId() {return spec_id_;}

  int getPeakId() {return peak_id_;}

  RmBreakTypePtr getBreakType();

  void setStrictTolerance(double tolerance){strict_tolerance_ = tolerance;}

  void setNStrictCRelacTolerance(double tolerance){
    n_strict_c_relax_tolerance_ = tolerance;}

  void setNRelaxCStrictTolerance(double tolerance){
    n_relax_c_strict_tolerance_ = tolerance;}

  void setPeakId(int peak_id) {peak_id_ = peak_id;}

  static bool cmpPosInc(const PrmPeakPtr &a, const PrmPeakPtr &b){
    return a->getPosition() < b->getPosition();
  }

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


} /* namespace prot */

#endif /* PRM_PEAK_HPP_ */
