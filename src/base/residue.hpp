/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_RESIDUE_HPP_
#define PROT_RESIDUE_HPP_

#include <string>
#include <memory>

#include "base/acid.hpp"
#include "base/ptm.hpp"

namespace prot {

class Residue {
 public:
	Residue(AcidPtr acid_ptr, PtmPtr ptm_ptr); 

  Residue(AcidPtrVec acid_list, PtmPtrVec ptm_list,
          std::string acid_name, std::string abbr_name);

	/** Get amino acid. */
	AcidPtr getAcidPtr() {return acid_ptr_; }
	
	/** Get residue mass. */
	double getMass() { return mass_; }

	/** Get post-translational modification. */
	PtmPtr getPtmPtr() { return ptm_ptr_; }

	/**
	 * Checks if the residue contains the same amino acid and ptm.
	 */
	bool isSame(AcidPtr acid_ptr, PtmPtr ptm_ptr) {
		return acid_ptr_.get() == acid_ptr.get() && ptm_ptr_.get() == ptm_ptr.get();
	}
	
	/** Get string representation */
  std::string toString(std::string delim_bgn, std::string delim_end);

 private:
	/** amino acid */
	AcidPtr acid_ptr_;
	/** post-translational modification */
	PtmPtr ptm_ptr_;
	/** residue mass */
	double mass_;
};

typedef std::shared_ptr<Residue> ResiduePtr;
typedef std::vector<ResiduePtr> ResiduePtrVec;

ResiduePtr getResiduePtrByAcid(ResiduePtrVec &residue_list,
                                  AcidPtr acid_ptr);
/**
 * Returns the first residue based on the acid and ptm. 
 */
ResiduePtr getResiduePtrByAcidPtm(ResiduePtrVec &residue_list,
                                  AcidPtr acid_ptr, PtmPtr ptm_ptr);

ResiduePtrVec getResiduePtrVecInstance(AcidPtrVec &acid_list, 
                                     PtmPtrVec &ptm_list,
                                     const char* file_name); 

ResiduePtr addResidue(ResiduePtrVec &residue_list, AcidPtr acid_ptr,
                      PtmPtr ptm_ptr);

ResiduePtrVec convertAcidToResidueSeq(ResiduePtrVec residue_list,
                                      AcidPtrVec acid_list);

}
#endif
