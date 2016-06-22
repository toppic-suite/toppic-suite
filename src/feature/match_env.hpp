#ifndef PROT_FEATURE_MATCH_ENVELOPE_HPP_
#define PROT_FEATURE_MATCH_ENVELOPE_HPP_

#include <memory>
#include <vector>

#include "feature/feature_mng.hpp" 
#include "feature/envelope.hpp" 
#include "feature/real_env.hpp" 

namespace prot {

class MatchEnv;

typedef std::shared_ptr<MatchEnv> MatchEnvPtr;

class MatchEnv {
 public:
  MatchEnv(int mass_group, EnvelopePtr theo_env_ptr, RealEnvPtr real_env_ptr);

  void compScr(FeatureMngPtr mng_ptr);

  static bool cmpScoreInc(const MatchEnvPtr &a, const MatchEnvPtr &b) { 
    return a->getScore() < b->getScore();
  }

  int getId() {return id_;}

  int getMassGroup() {return mass_group_;}

  RealEnvPtr getRealEnvPtr() {return real_env_ptr_;}

  EnvelopePtr getTheoEnv() {return theo_env_ptr_;}

  double getScore() {return score_;}

  void setId(int id) {id_ = id;}

  void setTheoEnvPtr(EnvelopePtr theo_env_ptr) {theo_env_ptr_ = theo_env_ptr;}

 private:
  int id_;
  // we divide envelopes into several groups based on monoisotopic  masses  
  int mass_group_;
  double score_;
  EnvelopePtr theo_env_ptr_;
  RealEnvPtr real_env_ptr_;

  double calcPeakScr(int id_x, double inte_sum, double tolerance);

  double calcShareInteAccu(int id_x, double inte_sum);

  double calcMzFactor(int id_x, double shift, double tolerance);

  double calcIntensityFactor(double theo_inte, double real_inte);

  double calcIntensityFactor(int id_x, double ratio);

  double findBestShift(FeatureMngPtr mng_ptr);

  double findBestRatio(FeatureMngPtr mng_ptr);

  double calcScrWithSftRatio(double shift, double ratio, double tolerance);
};

typedef std::vector<MatchEnvPtr> MatchEnvPtrVec;

}

#endif
