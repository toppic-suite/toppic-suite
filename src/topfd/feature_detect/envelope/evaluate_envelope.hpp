//
// Created by abbash on 9/5/22.
//

#ifndef TOPPIC_EVALUATE_ENVELOPE_HPP
#define TOPPIC_EVALUATE_ENVELOPE_HPP

#include "seed_envelope.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/env_set/env_set_util.hpp"
#include "topfd/feature_detect/util/utility_functions.hpp"
#include "topfd/feature_detect/score/env_util.hpp"

namespace toppic {
namespace evaluate_envelope {
    bool preprocess_env(PeakMatrix &peak_matrix, SeedEnvelope &seed_env, double mass_tol, double corr_tol, bool valid);
    bool evaluate_envelope(PeakMatrix &peak_matrix, SeedEnvelope &seed_envelope, double mass_tol, double corr_tol);
    bool test_charge_state(int charge, std::vector<double> &seed_envelope_inte);
    bool evaluate_envelope_pair(std::vector<double> &experimental_envelope_inte, std::vector<double> &theo_inte, double corr_tol);
}
}


#endif //TOPPIC_EVALUATE_ENVELOPE_HPP
