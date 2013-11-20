#ifndef PROTOMICS_RESIDUE_VEC_H_
#define PROTOMICS_RESIDUE_VEC_H_

#include <vector>

#include "residue.hpp"

namespace proteomics {

class ResidueVec {
 public:
  static ResidueVec(

	public ResList() throws Exception {
		acidList = AcidList.getCompleteAcidList();
		ptmList = PtmList.getCompletePtmList();
	}
	
	
	/**
	 * Adds a new residue into the list based on the name of the amino acid and
	 * the abbreviation name of the ptm.
	 */
	public Res add(String acidName, String ptmAbbrName) {
		Acid acid = acidList.getAcidByName(acidName);
		Ptm ptm = ptmList.getPtmByAbbrName(ptmAbbrName);
		return addRes(acid, ptm);
	}
	
	/** Adds a new residue */
	public Res addRes(Acid acid, Ptm ptm) {
		if (acid == null || ptm == null) {
			return null;
		}
		/* check if the residue is in the list */
		Res res = getResByAcidPtm(acid, ptm);
		if (res == null) {
			res = new Res(acid, ptm);
			add(res);
		}
		return res;
	}
	
	/**
	 * Returns the first residue based on the its acid. Returns null if the
	 * residue does not exist.
	 */
	public Res getResByAcid(Acid acid) {
		for (int i = 0; i < size(); i++) {
			if (get(i).getAcid() == acid) {
				return get(i);
			}
		}
		logger.debug("Acid not found " + acid.getOneLetter());
		return null;
	}

	/**
	 * Returns the first residue based on the one letter presentation of the
	 * amino acid.
	 */
	public Res getResByAcid(String acidOneLetter) {
		Acid acid = acidList.getAcidByOneLetter(acidOneLetter);
		return getResByAcid(acid);
	}

	/**
	 * Returns the first residue based on the acid and ptm. 
	 */
	public Res getResByAcidPtm(Acid acid, Ptm ptm) {
		for (int i = 0; i < size(); i++) {
			if (get(i).isSame(acid, ptm)) {
				return get(i);
			}
		}
		logger.debug("Residue not found " + acid.getOneLetter() + " " + ptm.getAbbrName());
		return null;
	}

	/** Get the first n residues in the list. */
	public Res[] getResidues(int n) {
        Res residues[] = new Res[n];
        for (int i = 0; i < n; i++) {
            residues[i] = get(i);
        }
        return residues;
	}
	
	public static ResList getCompleteList() {
		return completeList;
	}
	
	public static Res addToCompleteList(Acid acid, Ptm ptm) {
		Res res = completeList.getResByAcidPtm(acid, ptm);
		if (res == null) {
			res = new Res(acid, ptm);
			completeList.add(res);
		}
		return res;
	}
}

