//
// Created by abbash on 7/29/22.
//

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_SET_ENV_SET_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_SET_ENV_SET_HPP

#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/env/seed_env.hpp"
#include "topfd/ecscore/env/ms_map_env.hpp"
#include "topfd/ecscore/env_set/xic.hpp"

namespace toppic {

class EnvSet;
typedef std::shared_ptr<EnvSet> EnvSetPtr;

class EnvSet {
 public:
  EnvSet(const SeedEnvPtr seed_ptr, MsMapEnvPtrVec env_list,
         int start, int end, double noise_inte_level, double sn_ratio);

  int getStartSpecId() const { return start_spec_id_; }

  void setStartSpecId(int start_spec_id) { start_spec_id_ = start_spec_id; }

  void setSpecId(int start_spec_id, int end_spec_id);

  int getEndSpecId() const { return end_spec_id_; }

  void setEndSpecId(int end_spec_id) { end_spec_id_ = end_spec_id; }

  int getCharge() { return seed_ptr_->getCharge(); }

  double getMass() { return seed_ptr_->getMonoNeutralMass(); }

  int getBaseSpecId() { return seed_ptr_->getSpecId(); }

  MsMapEnvPtrVec getExpEnvList() { return exp_env_list_; }

  int countEnvNum();

  void setExpEnvList(MsMapEnvPtrVec exp_env_list) { exp_env_list_ = exp_env_list; }

  SeedEnvPtr getSeedPtr() { return seed_ptr_; }

  XicPtr getXicPtr() { return xic_ptr_; }

  double getXicSeedInte();

  void setXicPtr(XicPtr xic_ptr) { xic_ptr_ = xic_ptr; }

  std::vector<double> getXicTopThreeInteList() { return xic_ptr_->getTopThreeInteList(); }

  std::vector<double> getXicAllPeakInteList() { return xic_ptr_->getAllPeakInteList(); }

  std::vector<double> getSeedInteList() {return seed_ptr_->getInteList(); }

  std::vector<double> getSeedMzList() {return seed_ptr_->getMzList(); }

  std::vector<double> compExpInteSumList();

  void getWeightMzError(double &cur_weight, double &cur_weight_mz_error);

  std::vector<std::vector<double>> getScaledTheoIntes(double sn_ratio, 
                                                      double noise_inte);

  void refineFeatureBoundary();

  void removePeakData(MsMapPtr matrix_ptr);

  double compIntensity(double sn_ratio, double noise_inte);

  static bool cmpChargeInc(EnvSetPtr a, EnvSetPtr b) { return a->getCharge() < b->getCharge(); }

 private:
  void initMedianXic();

  void shortlistExpEnvs();

 private:
  SeedEnvPtr seed_ptr_;
  MsMapEnvPtrVec exp_env_list_;
  XicPtr xic_ptr_;
  int start_spec_id_;
  int end_spec_id_;
  double min_inte_;
};

typedef std::vector<EnvSetPtr> EnvSetPtrVec;

}

#endif //TOPPIC_ENV_SET_HPP
