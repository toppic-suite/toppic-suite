#include "base/logger.hpp"
#include "tdgf/comp_pvalue_array.hpp"

namespace prot {

CompPValueArray::CompPValueArray(ProteoformPtrVec &raw_forms, 
                                 ProteoformPtrVec &prot_mod_forms,
                                 ResFreqPtrVec &prot_n_term_residues,
                                 ResFreqPtrVec &pep_n_term_residues,
                                 ResFreqPtrVec &residues,
                                 TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  prot_comp_prob_ptr_ = CompProbValuePtr(
      new CompProbValue(mng_ptr_->double_to_int_constant_,
                        prot_n_term_residues, residues, 
                        mng_ptr_->unexpected_shift_num_ + 1, 
                        mng_ptr_->max_table_height_, mng_ptr_->max_sp_prec_mass_));
  pep_comp_prob_ptr_ = CompProbValuePtr(
      new CompProbValue(mng_ptr_->double_to_int_constant_,
                        pep_n_term_residues, residues, 
                        mng_ptr_->unexpected_shift_num_ + 1,
                        mng_ptr_->max_table_height_, mng_ptr_->max_sp_prec_mass_));
                        
  test_num_ptr_ = CountTestNumPtr(new CountTestNum(raw_forms, prot_mod_forms,
                                                   prot_n_term_residues, residues, mng_ptr));
}

/* set alignment */
ExtremeValuePtrVec CompPValueArray::compExtremeValue(PrmMsPtr ms_six, 
                                                     PrSMPtrVec &prsms, bool strict) {
  PrmPeakPtrVec prm_peaks = ms_six->getPeakPtrVec();
  std::vector<double> prot_probs; 
  compProbArray(prot_comp_prob_ptr_, prm_peaks, prsms, strict, prot_probs);
  std::vector<double> pep_probs;
  compProbArray(pep_comp_prob_ptr_, prm_peaks, prsms, strict, pep_probs);
  double prec_mass = ms_six->getHeaderPtr()->getPrecMonoMassMinusWater();
  double tolerance = ms_six->getHeaderPtr()->getErrorTolerance();
  ExtremeValuePtrVec ev_probs; 
  for (unsigned int i = 0; i < prsms.size(); i++) {
    int inter_shift_num = prsms[i]->getProteoformPtr()->getUnexpectedChangeNum();
    SemiAlignTypePtr t = prsms[i]->getProteoformPtr()->getSemiAlignType();
    double cand_num = test_num_ptr_->compCandNum(t, inter_shift_num, 
                                                 prec_mass, tolerance);
    if (cand_num == 0.0) {
      LOG_WARN("Zero candidate number");
      cand_num = std::numeric_limits<double>::infinity();
    }
    if (t == SemiAlignTypeFactory::getCompletePtr() 
        || t == SemiAlignTypeFactory::getPrefixPtr()) {
      ev_probs.push_back(ExtremeValuePtr(new ExtremeValue(prot_probs[i], cand_num, 1)));
      /*
      logger.debug("AlignType " + t + " shift "
                   + prsms[i].getAnnoProtein().getUnknownShiftNum()
                   + " thresh " + prsms[i].getMatchFragNum()
                   + " prot one prob " + protProbs[i] + " pep one prob "
                   + pepProbs[i] + " candiate num " + nCandidates
                   + " p value " + eVProbs[i].getPValue());
                   */
    } else {
      ev_probs.push_back(ExtremeValuePtr(new ExtremeValue(pep_probs[i], cand_num,1)));
      /*
      logger.debug("AlignType " + t + " shift "
                   + prsms[i].getAnnoProtein().getUnknownShiftNum()
                   + " thresh " + prsms[i].getMatchFragNum()
                   + " one prob " + pepProbs[i] + " candiate num "
                   + nCandidates + " p value " + eVProbs[i].getPValue());
                   */
    }
  }
  return ev_probs;
}

}

/*

    private ExtremeValueProb compExtremeValueProb(Ms<PrmPeak> msSix, PrSM prsm) {
        PrSM prsms[] = new PrSM[1];
        prsms[0] = prsm;
        ExtremeValueProb evProbs[] = compExtremeValue(msSix, prsms, true);
        return evProbs[0];
    }

    public void setPValue(Ms<DeconvPeak> deconvMs, PrSM prsm) throws Exception {
        double refinePrecMass = prsm.getAdjustedPrecMass();
        Ms<DeconvPeak> refineDeconvSp = DeconvSpFactory.getRefineMs(deconvMs,
                prsm.getCalibration(), refinePrecMass);

        logger.debug("recalibration " + prsm.getCalibration()
                + " original precursor "
                + deconvMs.getHeader().getPrecMonoMass() + " precursor "
                + refinePrecMass);
        // Delta = 0 is important 
        Ms<PrmPeak> prmSpSix = PrmSpFactory.getSpSix(refineDeconvSp, 0, mng.spPara);
        // AlignMs<PrmPeak> prmSpSix = AlignMsFactory.getSpSix(deconvMs, 0,
        // mng);
        ExtremeValueProb prob = compExtremeValueProb(prmSpSix, prsm);
        prsm.setProb(prob);
    }

    public void setPValueArray(Ms<PrmPeak> msSix, PrSM prsms[]) {
        ExtremeValueProb eVProbs[] = compExtremeValue(msSix, prsms, false);
        for (int i = 0; i < prsms.length; i++) {
            prsms[i].setProb(eVProbs[i]);
        }
    }
}
*/
