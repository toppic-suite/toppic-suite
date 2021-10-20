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

#ifndef TOPPIC_STAT_LOCAL_UTIL_HPP_
#define TOPPIC_STAT_LOCAL_UTIL_HPP_

#include "common/base/ptm.hpp"
#include "seq/mass_shift.hpp"
#include "prsm/prsm.hpp"
#include "stat/local/local_mng.hpp"

namespace toppic {

namespace local_util {

MassShiftPtrVec massShiftFilter(const MassShiftPtrVec & mass_shift_vec, AlterTypePtr type);

void scrFilter(std::vector<double> & scr, int & bgn, int & end, double & conf, double threshold);

PtmPtrVec getPtmPtrVecByMass(double mass, double err, const PtmPtrVec & ptm_vec);

PtmPairVec getPtmPairVecByMass(double mass, double err, const PtmPairVec & ptm_pair_vec);

void compMTable(const std::vector<double> &ori_theo_masses, double shift, 
    const std::vector<double> &exp_masses, 
    PeakTolerancePtr tole_ptr,  
    std::vector<int> &m_table); 

void addTwoVectors(std::vector<int> &row_1, std::vector<int> &row_2); 

std::vector<int> compPrefScore(std::vector<int> & row); 

std::vector<int> compSuffScore(std::vector<int> & row); 

void compOnePtmSTable(std::vector<int> &s_table, ProteoformPtr form_ptr, 
                      const ExtendMsPtrVec & extend_ms_ptr_vec,
                      PtmPtr ptm_ptr, LocalMngPtr mng_ptr); 

void compTwoPtmMTable(std::vector<std::vector<int>> &s_table, ProteoformPtr form_ptr, 
                      const ExtendMsPtrVec & extend_ms_ptr_vec,
                      PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2, LocalMngPtr mng_ptr); 

void normalize(std::vector<double> & scr);

int compMatchFragNum(ProteoformPtr proteoform_ptr, const ExtendMsPtrVec &ms_ptr_vec, double min_mass);

}  // namespace local_util

}  // namespace toppic

#endif
