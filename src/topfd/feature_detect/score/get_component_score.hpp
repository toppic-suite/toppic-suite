//
// Created by abbash on 8/26/22.
//

#ifndef TOPPIC_GET_COMPONENT_SCORE_HPP
#define TOPPIC_GET_COMPONENT_SCORE_HPP

#include "env_util.hpp"
#include "topfd/envcnn/env_cnn.hpp"
#include "topfd/feature_detect/env_collection/env_collection.hpp"

namespace toppic {
namespace component_score {
    double get_agg_odd_even_peak_ratio(EnvCollection& env_coll);
    double get_agg_env_corr(EnvCollection& env_coll);
    double get_3_scan_corr(EnvCollection& env_coll);
    double get_matched_peaks_percent(EnvCollection& env_coll, std::vector<std::vector<double>> theo_map);
    double get_num_theo_peaks(std::vector<std::vector<double>>& theo_map);
    double get_consecutive_peaks_percent(EnvCollection& env_coll);
    double get_rt_range(EnvCollection& env_coll);
    double get_charge_range(EnvCollection env_coll);
    double get_mz_errors(EnvCollection& env_coll);
    double count_max_consecutive_peak_num(std::vector<ExpEnvelope> & peaks);
}
}


#endif //TOPPIC_GET_COMPONENT_SCORE_HPP
