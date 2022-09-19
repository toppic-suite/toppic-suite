//
// Created by abbash on 8/26/22.
//

#include <cmath>
#include "utility_functions.hpp"

namespace toppic{
namespace utility_functions{
  double pearsonr(std::vector<double> & X, std::vector<double> & Y){
    int n = X.size();
    double sum_X = 0, sum_Y = 0, sum_XY = 0;
    double squareSum_X = 0, squareSum_Y = 0;
    for (int i = 0; i < n; i++){
      sum_X = sum_X + X[i];
      sum_Y = sum_Y + Y[i];
      sum_XY = sum_XY + X[i] * Y[i];
      squareSum_X = squareSum_X + X[i] * X[i];
      squareSum_Y = squareSum_Y + Y[i] * Y[i];
    }
    double corr = (double)(n * sum_XY - sum_X * sum_Y) / sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
    return corr;
  }

  std::vector<int> findLocalMinima(std::vector<double> & arr) {
    int n = arr.size();
    std::vector<int> minima;
    for (int i = 1; i < n - 1; i++) {
      if ((arr[i - 1] > arr[i]) and (arr[i] < arr[i + 1])) {
        if (i - 2 > 0)
          if (arr[i - 2] <= arr[i])
            continue;
        if (i + 2 < n)
          if (arr[i + 2] <= arr[i])
            continue;
        minima.push_back(i);
      }
    }
    return minima;
  }

  std::vector<int> findLocalMaxima(std::vector<double> & arr) {
    int n = arr.size();
    std::vector<int> maxima;
//    if (arr[0] > arr[1]) maxima.push_back(0);
    for (int i = 1; i < n - 1; i++)
      if ((arr[i - 1] < arr[i]) and (arr[i] > arr[i + 1]))
        maxima.push_back(i);
    if (arr[n - 1] > arr[n - 2]) maxima.push_back(n - 1);
    return maxima;
  }

  toppic::SeedEnvelope get_half_charge_env(SeedEnvelope &env, double even_odd_peak_ratios) {
    double mass = env.getMass();
    int charge = env.getCharge();
    double mz = env_utils::get_mz(mass, charge);
    std::vector<double> distribution = env.get_pos_list();

    if (even_odd_peak_ratios < 0)
      mz = mz + (distribution[1] - distribution[0]);
    int new_charge = int(charge / 2);
    if (new_charge == 0)
      new_charge = new_charge + 1;
    double new_mass = env_utils::get_mass(mz, new_charge);

    // get a reference distribution based on the base mass
    EnvelopePtr ref_env_ptr = EnvBase::getStaticEnvByMonoMass(new_mass);
    if (ref_env_ptr == nullptr) return SeedEnvelope();
    EnvelopePtr theo_env_ptr = ref_env_ptr->distrToTheoMono(mz, new_charge);
    std::vector<double> env_peaks_mz, env_peaks_inte;
    for (int i = 0; i < theo_env_ptr->getPeakNum(); i++) {
      env_peaks_mz.push_back(theo_env_ptr->getMz(i));
      env_peaks_inte.push_back(theo_env_ptr->getIntensity(i));
    }
    SeedEnvelope sp_peak = SeedEnvelope(env.getSpecId(), env.getEnvId(), env.getPos(), new_mass, env.getInte(), new_charge, env_peaks_mz, env_peaks_inte);
    return sp_peak;
  }

  toppic::SeedEnvelope test_half_charge_state(PeakMatrix &peak_matrix, SeedEnvelope &env, EnvSet &top_peak_env_set, double even_odd_peak_ratios, double mass_tole) {
    SeedEnvelope half_charge_env = get_half_charge_env(env, even_odd_peak_ratios);
    bool valid = false;
    valid = evaluate_envelope::preprocess_env(peak_matrix, half_charge_env, mass_tole, 0.5, valid);
    if (!valid) return SeedEnvelope();
    return half_charge_env;
  }



}
}