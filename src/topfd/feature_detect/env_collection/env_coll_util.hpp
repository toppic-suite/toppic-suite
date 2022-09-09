//
// Created by abbash on 8/24/22.
//

#ifndef TOPPIC_ENV_COLL_UTIL_HPP
#define TOPPIC_ENV_COLL_UTIL_HPP

#include "env_collection.hpp"
//#include "topfd/feature_detect/env_set/env_set_util.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "ms/env/env_base.hpp"

namespace toppic {
namespace env_coll_util {
  EnvCollection find_env_collection(PeakMatrix& peak_matrix, SeedEnvelope& seed_env, double mass_tole, int max_miss_env,
                                             int max_miss_charge, int max_miss_peak, int para_max_charge, double ratio_multi,
                                             double match_peak_tole, double snr);
  std::vector<EnvSet> get_charge_env_list(PeakMatrix& peak_matrix, SeedEnvelope& env, EnvSet top_peak_env_set, double para_max_charge, double mass_tole, int max_miss_peak, double max_miss_charge, double ratio_multi, double snr);
  void get_charge_range(const PeakMatrix& peak_matrix, SeedEnvelope seed_env, double mass_tole, int max_miss_num,
                        int para_max_charge, double match_peak_tole, int* return_min_charge, int* return_max_charge);
  SeedEnvelope select_best_seed_envelope(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tole, const EnvBase& env_base);

  void output_file();

}
}


#endif //TOPPIC_ENV_COLL_UTIL_HPP
