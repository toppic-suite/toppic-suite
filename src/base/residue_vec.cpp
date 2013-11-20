/**
 * Describes a list of residues.
 *
 * @author  Xiaowen Liu
 * @date    2009-8-26
 */

package com.jap.proteomics.base.residue;

import java.util.ArrayList;

import org.apache.log4j.Logger;

public class ResList extends ArrayList<Res> {

	private static Logger logger = Logger.getLogger(ResList.class);
	private static final long serialVersionUID = 1L;
	private AcidList acidList;
	private PtmList ptmList;
	
	/** Complete residue list, all other residue list are subset of this list */
	private static ResList completeList;
	static {
		try {
			completeList = ResListFactory.getDefinedList();
		}
		catch (Exception ex) {
			System.err.println("Error in initializing the complete ptm list");
			System.exit(1);
		}
	}

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

