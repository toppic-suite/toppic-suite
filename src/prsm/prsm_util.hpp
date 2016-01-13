#ifndef PROT_PRSM_PRSM_UTIL_HPP_
#define PROT_PRSM_PRSM_UTIL_HPP_

#include <memory>
#include <vector>
#include <string>

#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {

class PrsmUtil {
 public:
  static std::string getValueStr(std::string line);

  static std::string getXmlLine(const std::vector<std::string> &str_vec,
                                const std::string &property);

  static PrsmPtrVec selectSpeciesPrsms(const PrsmPtrVec &prsm_ptrs,int species_id);

  static std::vector<int> getSpeciesIds(const PrsmPtrVec &prsm_ptrs, std::string &seq_name);
  
  static int getProteinId(const PrsmPtrVec &prsm_ptrs, std::string &seq_name);

  static std::vector<int> getSpeciesIds(const PrsmPtrVec &prsm_ptrs);

  static void addSpectrumPtrsToPrsms(PrsmPtrVec &prsm_ptrs, PrsmParaPtr prsm_para_ptr);
};

}
#endif

