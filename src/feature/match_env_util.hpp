#ifndef PROT_FEATURE_MATCH_ENV_UTIL_HPP_
#define PROT_FEATURE_MATCH_ENV_UTIL_HPP_

#include "spec/peak.hpp"
#include "feature/match_env.hpp"

namespace prot {

class MatchEnvUtil {
 public:
  static MatchEnvPtrVec sortOnMz(MatchEnvPtrVec &ori_envs);

  static MatchEnvPtrVec sortOnMass(MatchEnvPtrVec &ori_envs);

  static std::vector<double> getMassList(MatchEnvPtrVec &envs);

  static std::vector<int> getChargeList(MatchEnvPtrVec &envs);

  static std::vector<double> getChargeOneMassList(MatchEnvPtrVec &envs);

  static std::vector<double> getIntensitySums(MatchEnvPtrVec &envs);

  static void assignIntensity(PeakPtrVec &ms, MatchEnvPtrVec &envs);

  static PeakPtrVec rmAnnoPeak(PeakPtrVec &ms, MatchEnvPtrVec &envs);

  static MatchEnvPtrVec addLowMassPeak(MatchEnvPtrVec &envs, std::vector<PeakPtr> &ms, 
                                       double tolerance);

  static MatchEnvPtr getNewMatchEnv(PeakPtrVec &ms, int idx, double tolerance);

  static MatchEnvPtrVec addMultipleMass(MatchEnvPtrVec &envs, MatchEnvPtr2D &candidates,
                                        double multi_min_mass, int multi_min_charge, double min_ratio);
  
};

}

#endif
