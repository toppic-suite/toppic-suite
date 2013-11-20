/*
 * author  Xiaowen Liu
 * date    2013-11-17
 */

#ifndef PROTOMICS_PTM_H_
#define PROTOMICS_PTM_H_

#include <string>
#include <vector>

#include "acid.hpp"

namespace proteomics {

class Ptm;
typedef std::shared_ptr<Ptm> PtmPtr;

class Ptm {
 public:
  Ptm(const std::string &abbr_name, 
      const std::vector<AcidPtr> &valid_acid_ptrs, 
      double mono_mass);

  /* Get amino acid composition. */
  std::string getAbbrName() { return abbr_name_;}

  /* Get  monoisotopic mass. */
  double getMonoMass() {return mono_mass_;}

  /* Get valid acid list. */
  std::vector<AcidPtr> getValidAcidPtrs() {return valid_acid_ptrs_;}

  /* Is it an empty PTM. */
  bool isEmpty();

  static PtmPtr getEmptyPtmPtr(std::vector<AcidPtr> &valid_acid_ptrs);

 private:
  /* Abbreviation name of a PTM */
  std::string abbr_name_;
  /* Valid acids of the PTM */
  std::vector<AcidPtr> valid_acid_ptrs_;
  /* monoisotopic mass */
  double mono_mass_;
};


}
#endif

