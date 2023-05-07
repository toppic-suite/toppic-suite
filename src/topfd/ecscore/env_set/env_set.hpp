//
// Created by abbash on 7/29/22.
//

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_SET_ENV_SET_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_SET_ENV_SET_HPP

#include "topfd/ecscore/envelope/seed_envelope.hpp"
#include "topfd/ecscore/envelope/exp_envelope.hpp"
#include "topfd/ecscore/env_set/xic.hpp"

namespace toppic {

class EnvSet;
typedef std::shared_ptr<EnvSet> EnvSetPtr;

class EnvSet {
 public:
  EnvSet(const SeedEnvelopePtr seed_ptr, ExpEnvelopePtrVec env_list, 
         int start, int end, double noise_inte_level, double sn_ratio);

  int getStartSpecId() const { return start_spec_id_; }

  void setStartSpecId(int start_spec_id) { start_spec_id_ = start_spec_id; }

  void setSpecId(int start_spec_id, int end_spec_id);

  int getEndSpecId() const { return end_spec_id_; }

  void setEndSpecId(int end_spec_id) { end_spec_id_ = end_spec_id; }

  int getCharge() { return seed_ptr_->getCharge(); }

  double getMass() { return seed_ptr_->getMass(); }

  int getBaseSpecId() { return xic_ptr_->getBaseSpecId(); }

  XicPtr getXicPtr() { return xic_ptr_; }

  void setXicPtr(XicPtr xic_ptr) { xic_ptr_ = xic_ptr; }

  static bool cmpCharge(EnvSetPtr a, EnvSetPtr b) { return a->getCharge() < b->getCharge(); }

 private:
  void initMedianXic(double noise_inte_level, double sn_ratio);

  /*
  void get_weight_mz_error(double *cur_weight, double *cur_weight_mz_error);

  std::vector<double> comp_exp_inte_sum_list();

  void refine_feature_boundary();

  std::vector<std::vector<double>> get_map(double snr, double noise_inte);

  double comp_intensity(double snr, double noise_inte);

  void remove_peak_data(PeakMatrix &peakMatrix);

  void shortlistExpEnvs();

  const toppic::SeedEnvelope getSeedEnv() const { return seed_env_; }

  void setSeedEnv(const toppic::SeedEnvelope &seedEnv) { seed_env_ = seedEnv; }

  const std::vector<ExpEnvelope> getExpEnvList() const { return exp_env_list_; }

  void setExpEnvList(const std::vector<ExpEnvelope> &expEnvList) { exp_env_list_ = expEnvList; }


  std::vector<double> getXicEnvIntes() { return xic_.getInteList(); }

  std::vector<double> getEnvIntes() { return xic_.getEnvInteList(); }


  std::vector<SimplePeak> get_peak_list() { return seed_env_.getPeakList(); }

  std::vector<double> get_theo_distribution_mz() { return seed_env_.get_pos_list(); }

  std::vector<double> get_theo_distribution_inte() { return seed_env_.get_inte_list(); }
  */



 private:
  SeedEnvelopePtr seed_ptr_;
  ExpEnvelopePtrVec exp_env_list_;
  XicPtr xic_ptr_;
  int start_spec_id_;
  int end_spec_id_;
};
}

#endif //TOPPIC_ENV_SET_HPP
