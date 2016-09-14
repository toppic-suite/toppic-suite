// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_ENVELOPE_HPP_
#define PROT_FEATURE_ENVELOPE_HPP_

#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include "spec/peak.hpp"

namespace prot {

class Envelope;

typedef std::shared_ptr<Envelope> EnvelopePtr;

class Envelope {
 public:

  Envelope() {}

  Envelope(const Envelope &env);

  Envelope(int num, std::vector<std::string> &line_List);

  Envelope(int refer_idx, int charge, double mono_mz,
           std::vector<double> &mzs, std::vector<double> &intensities);

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

  int getCharge() {return charge_;}

  int getLabel(int i) {return (int)std::round((mzs_[i] - mono_mz_) * charge_);}

  double getIntensity(int i) {return intensities_[i];}

  std::vector<double> getIntensities() {return intensities_;}

  double getMonoMass() {return Peak::compPeakMass(mono_mz_, charge_);}

  double getMonoMz() {return mono_mz_;}

  // get the m/z difference between mono_mz and reference peak 
  double getMonoReferDistance() {return mzs_[refer_idx_] - mono_mz_;}

  double getMz(int i) {return mzs_[i];}

  int getPeakNum() {return mzs_.size();}

  int getReferIdx() {return refer_idx_;}

  double getReferIntensity() {return intensities_[refer_idx_];}

  double getReferMz() {return mzs_[refer_idx_];}

  void setIntensity(int i, double intensity) {intensities_[i] = intensity;}

 protected:
  
  int refer_idx_;
  // Charge of the envolope 
  int charge_;
  // Theoretical m/z value of monoisotopic ion 
  double mono_mz_;
  // m/z list of all peaks in this envelope 
  std::vector<double> mzs_;
  // intensity list of all peaks in this envelope 
  std::vector<double> intensities_;
};

typedef std::vector<EnvelopePtr> EnvelopePtrVec;

}

#endif
