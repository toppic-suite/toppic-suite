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

  ExpEnvelopePtrVec getExpEnvList() { return exp_env_list_; }

  void setExpEnvList(ExpEnvelopePtrVec exp_env_list) { exp_env_list_ = exp_env_list; }

  SeedEnvelopePtr getSeedPtr() { return seed_ptr_; }

  XicPtr getXicPtr() { return xic_ptr_; }

  void setXicPtr(XicPtr xic_ptr) { xic_ptr_ = xic_ptr; }

  std::vector<double> getXicEnvInteList() { return xic_ptr_->getInteList(); }

  std::vector<double> getSeedInteList() {return seed_ptr_->getInteList(); }

  std::vector<double> getSeedMzList() {return seed_ptr_->getPosList(); }

  std::vector<double> compExpInteSumList();

  void getWeightMzError(double &cur_weight, double &cur_weight_mz_error);

  std::vector<std::vector<double>> getScaledTheoIntes(double sn_ratio, 
                                                      double noise_inte);

  void refineFeatureBoundary();

  double compIntensity(double sn_ratio, double noise_inte);

  static bool cmpCharge(EnvSetPtr a, EnvSetPtr b) { return a->getCharge() < b->getCharge(); }

 private:
  void initMedianXic(double noise_inte_level, double sn_ratio);

  void shortlistExpEnvs();

  /*



  std::vector<std::vector<double>> get_map(double snr, double noise_inte);


  void remove_peak_data(PeakMatrix &peakMatrix);



  void setSeedEnv(const toppic::SeedEnvelope &seedEnv) { seed_env_ = seedEnv; }



  std::vector<double> getXicEnvIntes() { return xic_.getInteList(); }

  std::vector<double> getEnvIntes() { return xic_.getEnvInteList(); }


  std::vector<SimplePeak> get_peak_list() { return seed_env_.getPeakList(); }

  std::vector<double> get_theo_distribution_mz() { return seed_env_.get_pos_list(); }

  */



 private:
  SeedEnvelopePtr seed_ptr_;
  ExpEnvelopePtrVec exp_env_list_;
  XicPtr xic_ptr_;
  int start_spec_id_;
  int end_spec_id_;
};

typedef std::vector<EnvSetPtr> EnvSetPtrVec;

}

#endif //TOPPIC_ENV_SET_HPP
