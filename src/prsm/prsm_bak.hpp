#ifndef PROT_PRSM_PRSM_FACTORY_HPP_
#define PROT_PRSM_PRSM_FACTORY_HPP_

#include <string>

#include "htslib/faidx.h"

#include "base/extreme_value.hpp"
#include "base/proteoform.hpp"
#include "spec/deconv_peak.hpp"
#include "spec/extend_peak.hpp"
#include "spec/sp_para.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {

class PrsmFactory {
 public:
  PrsmPtrVec readPrsm(const std::string &file_name, const ProteoformPtrVec &proteoforms);

void filterPrsms(const PrsmPtrVec &prsms, MsHeaderPtr header_ptr, PrsmPtrVec &sele_prsms); 

void addSpectrumPtrsToPrsms(PrsmPtrVec &prsms, PrsmParaPtr prsm_para_ptr);

PrsmPtrVec selectSpeciesPrsms(const PrsmPtrVec &prsm_ptrs,int species_id);

std::vector<int> getSpeciesIds(const PrsmPtrVec &prsm_ptrs,int seq_id);

std::vector<int> getSpeciesIds(const PrsmPtrVec &prsm_ptrs);

}
#endif

