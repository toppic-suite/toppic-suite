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


#ifndef PROT_PRSM_PRSM_UTIL_HPP_
#define PROT_PRSM_PRSM_UTIL_HPP_

#include <memory>
#include <vector>
#include <string>

#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {
namespace PrsmUtil {
std::string getValueStr(std::string line);

std::string getXmlLine(const std::vector<std::string> &str_vec,
                       const std::string &property);

std::vector<std::string> getXmlLineVec(const std::vector<std::string> &str_vec,
                                       const std::string &property);

PrsmPtrVec selectSpeciesPrsms(const PrsmPtrVec &prsm_ptrs, int species_id);

std::vector<int> getSpeciesIds(const PrsmPtrVec &prsm_ptrs, std::string &seq_name);

int getProteinId(const PrsmPtrVec &prsm_ptrs, std::string &seq_name);

std::vector<int> getSpeciesIds(const PrsmPtrVec &prsm_ptrs);

void addSpectrumPtrsToPrsms(PrsmPtrVec &prsm_ptrs, PrsmParaPtr prsm_para_ptr);
}  // namesace PrsmUtil
}  // namespace prot
#endif

