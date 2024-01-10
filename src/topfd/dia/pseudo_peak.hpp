//
// Created by abbash on 12/30/23.
//

#ifndef TOPPIC_TOPFD_DIA_PSEUDO_PEAK_HPP
#define TOPPIC_TOPFD_DIA_PSEUDO_PEAK_HPP

namespace toppic {

class PseudoPeak {
  public:

  PseudoPeak(double mass, double mono_mz, int charge, double intensity, double score, double corr, double coverage,
             double apex_diff, double apex_diff_scan, double rt_low, double rt_high); 

  double getMass() const { return mass_; }

  double getMonoMz() const { return mono_mz_; }

  int getCharge() const { return charge_; }

  double getIntensity() const { return intensity_; }

  double getScore() const { return score_; }

  double getCorr() const { return corr_; }

  double getCoverage() const { return coverage_; }

  double getApexDiff() const { return apex_diff_; }

  double getApexDiffScan() const { return apex_diff_scan_; }

  double getRtLow() const { return rt_low_; }

  double getRtHigh() const { return rt_high_; }

  private:
  double mass_;
  double mono_mz_;
  int charge_;
  double intensity_;
  double score_;
  double corr_;
  double coverage_;
  double apex_diff_;
  double apex_diff_scan_;
  double rt_low_;
  double rt_high_;
};
}

#endif 
