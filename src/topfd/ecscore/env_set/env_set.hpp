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

  EnvSet(const SeedEnvPtr seed_ptr, MsMapEnvPtrVec env_list,
         int start, int end, double min_inte);

  int getStartSpecId() const { return start_spec_id_; }

  void setStartSpecId(int start_spec_id) { start_spec_id_ = start_spec_id; }

  int getEndSpecId() const { return end_spec_id_; }

  void setEndSpecId(int end_spec_id) { end_spec_id_ = end_spec_id; }

  int getCharge() { return seed_ptr_->getCharge(); }

  double getMonoMass() { return seed_ptr_->getMonoNeutralMass(); }

  int getSeedSpecId() { return seed_ptr_->getSpecId(); }

  MsMapEnvPtrVec getMsMapEnvList() { return ms_map_env_list_; }

  int countEnvNum();

  void setMsMapEnvList(MsMapEnvPtrVec ms_map_env_list) { ms_map_env_list_ = ms_map_env_list; }

  SeedEnvPtr getSeedPtr() { return seed_ptr_; }

  void setSeedPtr(SeedEnvPtr seed_ptr) {seed_ptr_ = seed_ptr;}

  XicPtr getXicPtr() { return xic_ptr_; }

  void setXicPtr(XicPtr xic_ptr) { xic_ptr_ = xic_ptr; }

  // get the all peak intensity in xic for the seed spectrum
  double getXicSeedAllPeakInte();

  // compute aggregate envelope peak intensities
  std::vector<double> compAggrEnvInteList();

  std::vector<double> compAggrEnvMzList(); 

  // seed peak intensity list x spectrum intensity ratio list
  std::vector<std::vector<double>> getScaledTheoIntes(int min_inte);

  double getInte() {return xic_ptr_->getAllPeakInteSum(); }

  void removePeakData(MsMapPtr matrix_ptr);

  std::pair<double, double> getMzErrorAndWeight();

  void refineXicBoundary();

  bool containTwoValidEnvs(int min_match_peak_num); 

  bool containValidNeighborEnvsForSeed(int min_match_peak_num);

  bool containThreeValidOutOfFiveEnvs(int min_match_peak_num); 

  void mergeEnvSet(EnvSetPtr new_set_ptr); 

  static bool cmpChargeInc(EnvSetPtr a, EnvSetPtr b) { return a->getCharge() < b->getCharge(); }

 private:
  void initMedianXic();

  void shortlistExpEnvs();

 private:
  SeedEnvPtr seed_ptr_;
  MsMapEnvPtrVec ms_map_env_list_;
  XicPtr xic_ptr_;
  int start_spec_id_;
  int end_spec_id_;
  double min_inte_;
};

typedef std::vector<EnvSetPtr> EnvSetPtrVec;

}

#endif //TOPPIC_ENV_SET_HPP
