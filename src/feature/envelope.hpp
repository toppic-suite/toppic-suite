#ifndef PROT_FEATURE_ENVELOPE_HPP_
#define PROT_FEATURE_ENVELOPE_HPP_

#include <memory>
#include <vector>
#include <string>

#include "spec/peak.hpp"

namespace prot {

class Envelope;

typedef std::shared_ptr<Envelope> EnvelopePtr;

class Envelope {
 public:

  Envelope() {}

  Envelope(const Envelope &env);

  Envelope(int num, std::string str);

  Envelope(int refer_idx, int charge, double mono_mz,
           std::vector<double> &mzs, std::vector<double> &intensities);

  EnvelopePtr convertToTheo(double mass_diff, int new_charge);

  EnvelopePtr distrToTheoBase(double new_base_mz, int new_charge);

  EnvelopePtr distrToTheoMono(double new_mono_mz, int new_charge);

  void changeIntensity(double ratio);

  void changeToAbsInte(double absolute_intensity);

  void changeMz(double shift);

  EnvelopePtr getSubEnv(int n_back, int n_forw);

  EnvelopePtr addZero(int num);

  EnvelopePtr getSubEnv(double percent_bound, double absolute_min_inte,
                        int max_back_peak_num, int max_forw_peak_num);

  std::vector<int> calcBound(double percent_bound, double absolute_min_inte,
                           int max_back_peak_num, int max_forw_peak_num);

  void shift(int shift);

  double compIntensitySum();

  double getAvgMz();

  double getAvgMass();

  int getCharge() {return charge_;}

  int getLabel(int i) {return (int)std::round((mzs_[i] - mono_mz_) * charge_);}

  double getIntensity(int i) {return intensities_[i];}

  std::vector<double> getIntensities() {return intensities_;}

  double getMonoMass() {return Peak::compPeakMass(mono_mz_, charge_);}

  double getMonoMz() {return mono_mz_;}

  // get the m/z difference between mono_mz and reference peak 
  double getMonoReferDistance() {return mzs_[refer_idx_] - mono_mz_;}

  double getMz(int i) {return mzs_[i];}

  int getPeakNum() {return mzs_.size();}

  int getReferIdx() {return refer_idx_;}

  double getReferIntensity() {return intensities_[refer_idx_];}

  double getReferMz() {return mzs_[refer_idx_];}

  void setIntensity(int i, double intensity) {intensities_[i] = intensity;}

 private:
  
  int refer_idx_;
  // Charge of the envolope 
  int charge_;
  // Theoretical m/z value of monoisotopic ion 
  double mono_mz_;
  // m/z list of all peaks in this envelope 
  std::vector<double> mzs_;
  // intensity list of all peaks in this envelope 
  std::vector<double> intensities_;

  std::vector<Peak> peak_list_;
};

typedef std::vector<EnvelopePtr> EnvelopePtrVec;

}

#endif
