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


#ifndef PROT_LOCAL_UTIL_HPP_
#define PROT_LOCAL_UTIL_HPP_

#include <vector>

#include "common/base/amino_acid_base.hpp"
#include "common/base/ptm.hpp"
#include "seq/mass_shift.hpp"
#include "common/base/ptm.hpp"

#include "prsm/prsm.hpp"
#include "prsm/peak_ion_pair_util.hpp"

#include "local/local_mng.hpp"

namespace toppic {

namespace local_util {

void scrFilter(std::vector<double> & scr, int & bgn, int & end, double & conf, double threshold);

PtmPtrVec getPtmPtrVecByMass(double mass, double err, const PtmPtrVec & ptm_vec);

PtmPairVec getPtmPairVecByMass(double mass1, double mass2, double err, const PtmPairVec & ptm_pair_vec);

void compSupPeakNum(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                    MassShiftPtr mass_shift, double min_mass, int & left, int & right);

void ptmMassAdjust(double & mass1, double & mass2, PtmPtr p1, PtmPtr p2);

void fillTableB(std::vector<std::vector<double>> & b_table, double mass1, double mass2);

void compNumMatch(const std::vector<double> & b, std::vector<int> & s,
                  const ExtendMsPtr & extend_ms_ptr, double prec_mass);

void fillTableS(std::vector<std::vector<double> > & b_table,
                std::vector<std::vector<int> > & s_table,
                ExtendMsPtr extend_ms_ptr, double prec_mass);

std::vector<double> geneNTheoMass(ProteoformPtr proteoform, ExtendMsPtr extend_ms_ptr_vec,
                                  double min_mass);

MassShiftPtrVec massShiftFilter(const MassShiftPtrVec & mass_shift_vec, MassShiftTypePtr type);

MassShiftPtrVec copyMassShiftVec(const MassShiftPtrVec & mass_shift_vec);

double compMassShift(const MassShiftPtrVec & mass_shift_vec);

MassShiftPtr geneMassShift(MassShiftPtr shift, double mass, MassShiftTypePtr type);

void normalize(std::vector<double> & scr);

int compMatchFragNum(ProteoformPtr proteoform_ptr, const ExtendMsPtrVec &ms_ptr_vec, double min_mass);

}  // namespace local_util

}  // namespace toppic

#endif
