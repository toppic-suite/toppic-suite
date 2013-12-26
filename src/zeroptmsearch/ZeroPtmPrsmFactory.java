package edu.ucsd.proteomics.search.zeroptmsearch;

import com.jap.proteomics.base.residue.EnumProtTermMod;
import com.jap.proteomics.base.theosp.BpSpec;
import com.jap.proteomics.base.theosp.NModBpSpec;

import edu.ucsd.proteomics.prsm.anno.AnnoProtein;
import edu.ucsd.proteomics.prsm.anno.EnumMassShiftType;
import edu.ucsd.proteomics.prsm.anno.MassShift;
import edu.ucsd.proteomics.prsm.base.SemiAlignType;

public class ZeroPtmPrsmFactory {
    public static AnnoProtein getAnnoProtein(NModBpSpec modSeq, int seqBgn,
            int seqEnd, SemiAlignType searchType) {
        EnumProtTermMod nMod = modSeq.getNMod();
        BpSpec unModSeq = modSeq.getUnModBpSpec();
        if (nMod == EnumProtTermMod.NONE) {
            return new AnnoProtein(unModSeq, seqBgn, seqEnd, new MassShift[0],
                    false, 0);
        } else if (nMod == EnumProtTermMod.ACETYLATION) {
            MassShift massShifts[] = new MassShift[1];
            massShifts[0] = new MassShift(-1, 1, EnumProtTermMod.ACETYLATION
                    .getPepTermShift(), EnumMassShiftType.EXPECTED_SITE);
            return new AnnoProtein(unModSeq, seqBgn, seqEnd, massShifts, true, 1);
        } else if (nMod == EnumProtTermMod.NME) {
            return new AnnoProtein(unModSeq, seqBgn + 1, seqEnd + 1, 
                    new MassShift[0], false, 0);
        } else if (nMod == EnumProtTermMod.NME_ACETYLATION) {
            MassShift massShifts[] = new MassShift[1];
            massShifts[0] = new MassShift(0, 2, EnumProtTermMod.NME_ACETYLATION
                    .getPepTermShift(), EnumMassShiftType.EXPECTED_SITE);
            return new AnnoProtein(unModSeq, seqBgn + 1, seqEnd + 1, 
                    massShifts, true, 1);
        }
        return null;
    }
}
