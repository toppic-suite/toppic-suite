//
// Created by abbash on 7/29/22.
//

#ifndef TOPPIC_ENV_SET_HPP
#define TOPPIC_ENV_SET_HPP

#include <vector>
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/spectrum/spectrum.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "exp_envelope.hpp"
#include "xic.hpp"

namespace toppic {
    class EnvSet {
    public:
        EnvSet();
        EnvSet(const SeedEnvelope& envelope, std::vector<ExpEnvelope> env_list, int start, int end);
        EnvSet(const EnvSet & es);

        void get_coordinates(spec_list spectra_list, std::vector<double> x, std::vector<double> y, std::vector<double> z);
        double get_median_ratio(ExpEnvelope env);
        Xic init_median_xic();
        double get_seed_inte_ratio();
        void remove_non_consecutive_peaks(int i, int max_miss_peak);
        void remove_all_non_consecutive_peaks(int max_miss_peak);
        void simple_remove_matrix_peaks(PeakMatrix peak_matrix);
        void remove_matrix_peaks(PeakMatrix peak_matrix);
        double comp_intensity(){ return seed_env_.get_inte_sum() * xic_.get_inte_list_sum();}
        void get_weight_mz_error(double* cur_weight, double* cur_weight_mz_error);
        std::vector<double> comp_exp_inte_sum_list();

        const toppic::SeedEnvelope getSeedEnv() const { return seed_env_; }
        void setSeedEnv(const toppic::SeedEnvelope &seedEnv) { seed_env_ = seedEnv; }

        const std::vector<ExpEnvelope> getExpEnvList() const { return exp_env_list_; }
        void setExpEnvList(const std::vector<ExpEnvelope> &expEnvList) { exp_env_list_ = expEnvList; }

        int getStartSpecId() const { return start_spec_id_; }
        void setStartSpecId(int startSpecId) { start_spec_id_ = startSpecId; }

        int getEndSpecId() const { return end_spec_id_; }
        void setEndSpecId(int endSpecId) { end_spec_id_ = endSpecId; }

        const Xic getXic() const { return xic_; }
        void setXic(const Xic &xic) { xic_ = xic; }

        std::vector<SimplePeak> get_peak_list() { return seed_env_.getPeakList(); }
        std::vector<double> get_theo_distribution_mz() {  return seed_env_.get_pos_list(); }
        std::vector<double> get_theo_distribution_inte() { return seed_env_.get_inte_list(); }

        bool isEmpty(){
          if (seed_env_.isEmpty() && exp_env_list_.size() == 0 && start_spec_id_ == -1 && end_spec_id_ == -1 && xic_.isEmpty())
            return true;
          return false;
        }

    private:
        SeedEnvelope seed_env_;
        std::vector<ExpEnvelope> exp_env_list_;
        int start_spec_id_;
        int end_spec_id_;
        Xic xic_;

    };
}

#endif //TOPPIC_ENV_SET_HPP
