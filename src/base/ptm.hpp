/*
 * author  Xiaowen Liu
 * date    2013-11-17
 */

#ifndef PROT_PTM_H_
#define PROT_PTM_H_

#include <string>
#include <vector>
#include <memory>

namespace prot {

class Ptm;
typedef std::shared_ptr<Ptm> PtmPtr;
typedef std::vector<PtmPtr> PtmPtrVec;

class Ptm {
 public:
  Ptm(const std::string &abbr_name, 
      double mono_mass);

  /* Get amino acid composition. */
  std::string getAbbrName() { return abbr_name_;}

  /* Get  monoisotopic mass. */
  double getMonoMass() {return mono_mass_;}

  /* Is it an empty PTM. */
  bool isEmpty();

  static PtmPtr getEmptyPtmPtr();

 private:
  /* Abbreviation name of a PTM */
  std::string abbr_name_;
  /* monoisotopic mass */
  double mono_mass_;
};

PtmPtrVec getPtmPtrVecInstance(const char* file_name);
/**
 * Returns a PTM based on the abbreviation name. Returns null if the
 * abbreviation name does not exist.
 */
PtmPtr getPtmPtrByAbbrName(PtmPtrVec &ptm_ptr_vec, 
                           const std::string &abbr_name);

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool containAbbrsName(PtmPtrVec &ptm_ptr_vec, const std::string &abbr_name);

PtmPtr findEmptyPtmPtr(PtmPtrVec &ptm_ptr_vec);

PtmPtr addPtm(PtmPtrVec &ptm_ptr_vec, std::string abbr_name,
              double mono_mass);

}
#endif

