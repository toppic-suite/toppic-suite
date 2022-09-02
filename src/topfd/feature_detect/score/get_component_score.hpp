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
    double get_agg_odd_even_peak_ratio(const EnvCollection& env_coll);
    double get_agg_env_corr(const EnvCollection& env_coll);
    double get_3_scan_corr(const EnvCollection& env_coll);
    double get_matched_peaks_percent(const EnvCollection& env_coll);
    double get_consecutive_peaks_percent(const EnvCollection& env_coll);
    double get_rt_range(const EnvCollection& env_coll);
    double get_charge_range(EnvCollection env_coll);
}
}


#endif //TOPPIC_GET_COMPONENT_SCORE_HPP
