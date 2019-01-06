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

#ifndef TOPPIC_PRSM_PRSM_UTIL_HPP_
#define TOPPIC_PRSM_PRSM_UTIL_HPP_

#include "prsm/prsm.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_para.hpp"

namespace toppic {

namespace prsm_util {

std::string getValueStr(std::string line);

std::string getXmlLine(const std::vector<std::string> &str_vec,
                       const std::string &property);

std::vector<std::string> getXmlLineVec(const std::vector<std::string> &str_vec,
                                       const std::string &property);

PrsmPtrVec selectClusterPrsms(const PrsmPtrVec &prsm_ptrs, int cluster_id);

std::vector<int> getProteoClusterIds(const PrsmPtrVec &prsm_ptrs, const std::string & seq_name);

int getProteinId(const PrsmPtrVec &prsm_ptrs, std::string &seq_name);

std::vector<int> getClusterIds(const PrsmPtrVec &prsm_ptrs);

void addSpectrumPtrsToPrsms(PrsmPtrVec &prsm_ptrs, PrsmParaPtr prsm_para_ptr);

void addFeatureIDToPrsms(PrsmStrPtrVec &prsm_ptrs, const std::string & feature_file_name);

void mergePrsmFiles(const std::vector<std::string> & prsm_file_lst, int N,
                    const std::string & output_file);

}  // namespace prsm_util

}  // namespace toppic
#endif

