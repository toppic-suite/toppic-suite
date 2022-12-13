//
// Created by abbash on 7/29/22.
//

#ifndef TOPPIC_ENV_SET_HPP
#define TOPPIC_ENV_SET_HPP

#include <vector>
#include <numeric>
#include <algorithm>
#include "topfd/feature_detect/env_set/exp_envelope.hpp"
#include "topfd/feature_detect/env_set/xic.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/spectrum/spectrum.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"


namespace toppic {
  class EnvSet;

  typedef std::shared_ptr<EnvSet> EnvSetPtr;

  class EnvSet {
  public:
    EnvSet();

    EnvSet(const SeedEnvelope &envelope, std::vector<ExpEnvelope> &env_list, int start, int end,
           double noise_inte_level, double snr);

    EnvSet(const EnvSet &es);

    Xic init_median_xic(double noise_inte_level, double snr);

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

    int getStartSpecId() const { return start_spec_id_; }

    void setStartSpecId(int startSpecId) { start_spec_id_ = startSpecId; }

    void setSpecId(int startSpecId, int endSpecId);

    int getEndSpecId() const { return end_spec_id_; }

    void setEndSpecId(int endSpecId) { end_spec_id_ = endSpecId; }

    int getBaseSpecId() { return xic_.getBaseSpecId(); }

    const Xic getXic() const { return xic_; }

    void setXic(const Xic &xic) { xic_ = xic; }

    std::vector<double> getXicEnvIntes() { return xic_.getInteList(); }

    std::vector<double> getEnvIntes() { return xic_.getEnvInteList(); }

    int getCharge() { return seed_env_.getCharge(); }

    double getMass() { return seed_env_.getMass(); }

    std::vector<SimplePeak> get_peak_list() { return seed_env_.getPeakList(); }

    std::vector<double> get_theo_distribution_mz() { return seed_env_.get_pos_list(); }

    std::vector<double> get_theo_distribution_inte() { return seed_env_.get_inte_list(); }

    bool isEmpty();

    static bool cmpCharge(EnvSet a, EnvSet b) { return a.getCharge() < b.getCharge(); }

  private:
    SeedEnvelope seed_env_;
    std::vector<ExpEnvelope> exp_env_list_;
    int start_spec_id_;
    int end_spec_id_;
    Xic xic_;
  };
}

#endif //TOPPIC_ENV_SET_HPP
