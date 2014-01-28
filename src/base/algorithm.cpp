#include <cmath>
#include <limits>

#include "base/algorithm.hpp"

namespace prot {

/* if we need to increase i, return true, otherwise, return false. */
bool increaseIJ(unsigned int i, unsigned int j, double deviation, 
                double tolerance, std::vector<double> ms_masses, 
                std::vector<double> theo_masses) {
        /*
         * we assume that each real peak is matched to at most one theoretical
         * peak, so we do not check i and j+1
         */
        if (deviation <= 0) {
            return true;
        }
        /* severl real peak can be matched to the same theoretical peak */
        if (i >= ms_masses[ms_masses.size() - 1] ) {
            return false;
        }

        double next_pos = ms_masses[i+1];
    bool j_is_near
        = std::abs(next_pos - theo_masses[j]) < std::abs(next_pos - theo_masses[j+1]);
    if (std::abs(next_pos - theo_masses[j]) <= tolerance  
        && (j == theo_masses.size() - 1 || j_is_near)) { 
            return true;
        } else {
            return false;
        }
}

/* compute deviation for each peak */
void  compMsMassPpos(std::vector<double> &ms_masses, 
                     std::vector<double> &theo_masses, 
                     double ppo, std::vector<double> &result_ppos) {
  // extendMsThree do not have 0 and precursor mass 
  std::vector<double> min_distances;
  for (unsigned p = 0; p < ms_masses.size(); p++) {
    min_distances.push_back(std::numeric_limits<double>::infinity());
  }
  unsigned int i = 0;
  unsigned int j = 0;
  while (i < ms_masses.size() && j < theo_masses.size()) {
    double d = ms_masses[i] - theo_masses[j];
    if (std::abs(d) <= std::abs(min_distances[i])) {
      min_distances[i] = d;
    }
    double tolerance = ms_masses[i] * ppo;
    if (increaseIJ(i, j, d,  tolerance, ms_masses, theo_masses)) {
      i++;
    } else {
      j++;
    }
  }
  // change distance to ppo
  for (i = 0; i < ms_masses.size(); i++) {
    if (ms_masses[i] > 0) {
      result_ppos.push_back (min_distances[i]/ ms_masses[i]); 
    }
    else {
      result_ppos.push_back(std::numeric_limits<double>::infinity());
    }
  }
}

void compTheoMassPpo(std::vector<double> &ms_masses, 
                     std::vector<double> &theo_masses,
                     double ppo, std::vector<double> &result_ppos) {

  std::vector<double> min_distances;
  for (unsigned p = 0; p < theo_masses.size(); p++) {
    min_distances.push_back(std::numeric_limits<double>::infinity());
  }
  /* extendMsThree do not have 0 and precursor mass */
  unsigned int i = 0;
  unsigned int j = 0;
  while (i < ms_masses.size() && j < theo_masses.size()) {
    double d = ms_masses[i] - theo_masses[j];
    if (std::abs(d) <= std::abs(min_distances[j])) {
      min_distances[j] = d;
    }
    double tolerance = ms_masses[i] * ppo;
    if (increaseIJ(i, j, d, tolerance, ms_masses, theo_masses)) {
      i++;
    } else {
      j++;
    }
  }
  // change distance to ppo
  for (i = 0; i < theo_masses.size(); i++) {
    if (theo_masses[i] > 0) {
      result_ppos.push_back (min_distances[i]/ theo_masses[i]); 
    }
    else {
      result_ppos.push_back(std::numeric_limits<double>::infinity());
    }
  }
}

/** compute unique score */
double compUniqueScore (std::vector<double> &ms_masses, 
                        std::vector<double> &theo_masses, 
                        double ppo) {
  std::vector<double> theo_mass_ppos;
  compTheoMassPpo(ms_masses, theo_masses, ppo, theo_mass_ppos);
  double score = 0;
  for (unsigned i = 0; i < theo_mass_ppos.size(); i++) {
    if (std::abs(theo_mass_ppos[i]) <= ppo) {
      score += 1.0;
    }
  }
  return score;
}

}
