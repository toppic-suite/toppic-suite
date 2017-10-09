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

#include <set>
#include <vector>
#include <string>
#include <algorithm>

#include "base/logger.hpp"
#include "spec/extend_ms_factory.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/prsm_util.hpp"

namespace prot {

std::string PrsmUtil::getValueStr(std::string line) {
  int start = line.find(">");
  int end = line.find("<", start);
  std::string num_str = line.substr(start + 1, end - start - 1);
  // LOG_DEBUG(line << "  Num str: " << num_str);
  return num_str;
}

std::string PrsmUtil::getXmlLine(const std::vector<std::string> &str_vec,
                                 const std::string &property) {
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      return str_vec[i];
    }
  }
  return "";
}

std::vector<std::string> PrsmUtil::getXmlLineVec(const std::vector<std::string> &str_vec,
                                                 const std::string &property) {
  std::vector<std::string> res;
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      res.push_back(str_vec[i]);
    }
  }
  return res;
}


PrsmPtrVec PrsmUtil::selectSpeciesPrsms(const PrsmPtrVec &prsm_ptrs, int species_id) {
  PrsmPtrVec select_prsm_ptrs;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (species_id == prsm_ptrs[i]->getProteoformPtr()->getSpeciesId()) {
      select_prsm_ptrs.push_back(prsm_ptrs[i]);
    }
  }
  return select_prsm_ptrs;
}

std::vector<int> PrsmUtil::getSpeciesIds(const PrsmPtrVec &prsm_ptrs, std::string &seq_name) {
  std::set<int> species_id_set;
  std::vector<int> species_ids;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getProteoformPtr()->getSeqName() == seq_name)
      species_id_set.insert(prsm_ptrs[i]->getProteoformPtr()->getSpeciesId());
  }
  std::copy(species_id_set.begin(), species_id_set.end(), std::back_inserter(species_ids));
  std::sort(species_ids.begin(), species_ids.end());
  return species_ids;
}

int PrsmUtil::getProteinId(const PrsmPtrVec &prsm_ptrs, std::string &seq_name) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getProteoformPtr()->getSeqName() == seq_name) {
      return prsm_ptrs[i]->getProteoformPtr()->getProtId();
    }
  }
  return -1;
}

std::vector<int> PrsmUtil::getSpeciesIds(const PrsmPtrVec &prsm_ptrs) {
  std::set<int> species_id_set;
  std::vector<int> species_ids;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    species_id_set.insert(prsm_ptrs[i]->getProteoformPtr()->getSpeciesId());
  }
  std::copy(species_id_set.begin(), species_id_set.end(), std::back_inserter(species_ids));
  std::sort(species_ids.begin(), species_ids.end());
  return species_ids;
}

bool isMatchMs(PrsmPtr prsm_ptr, MsHeaderPtr header_ptr) {
  int id = header_ptr->getId();
  std::string scan = header_ptr->getScansString();
  int prec_id = header_ptr->getPrecId();
  if (prsm_ptr->getSpectrumId() == id && prsm_ptr->getPrecursorId() == prec_id) {
    if (prsm_ptr->getSpectrumScan() != scan) {
      LOG_ERROR("Error in Prsm.");
    }
    return true;
  } else {
    LOG_TRACE("prsm spectrum id " << prec_id << " ms spectrum id " << id);
    return false;
  }
}

void PrsmUtil::addSpectrumPtrsToPrsms(PrsmPtrVec &prsm_ptrs, PrsmParaPtr prsm_para_ptr) {
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(),
                       prsm_para_ptr->getGroupSpecNum(),
                       prsm_para_ptr->getSpParaPtr()->getActivationPtr(),
                       prsm_para_ptr->getSpParaPtr()->getSkipList());
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  SpectrumSetPtr spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)[0];
  // use prsm order information (ordered by spectrum id then prec id)
  int start_prsm = 0;
  while (spec_set_ptr != nullptr) {
    if (spec_set_ptr->isValid()) {
      DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
      int spectrum_id = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getId();
      int prec_id = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecId();
      LOG_TRACE("spectrum id " << spectrum_id);
      for (size_t i = start_prsm; i < prsm_ptrs.size(); i++) {
        if (isMatchMs(prsm_ptrs[i], deconv_ms_ptr_vec[0]->getMsHeaderPtr())) {
          prsm_ptrs[i]->setDeconvMsPtrVec(deconv_ms_ptr_vec);
          double new_prec_mass = prsm_ptrs[i]->getAdjustedPrecMass();
          prsm_ptrs[i]->setRefineMsVec(ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec,
                                                                          prsm_para_ptr->getSpParaPtr(),
                                                                          new_prec_mass));
        }
        if ((spectrum_id == prsm_ptrs[i]->getSpectrumId() && prec_id < prsm_ptrs[i]->getPrecursorId()) ||
            spectrum_id < prsm_ptrs[i]->getSpectrumId()) {
          start_prsm = i;
          break;
        }
      }
    }
    spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)[0];
  }
  reader.close();
}

}  // namespace prot
