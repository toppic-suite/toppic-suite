#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "feature/env_dist.hpp" 

namespace prot {

EnvDist::EnvDist(std::string file_name, int entry_num, double mass_interval):
    entry_num_(entry_num),
    mass_interval_(mass_interval) {
      LOG_DEBUG("Distribution envelope file name " << file_name);

      std::ifstream input;
      input.open(file_name.c_str(), std::ios::in);
      for (int i = 0; i < entry_num_; i++) {
        int peak_num = 0;
        std::string line;
        std::vector<std::string> line_list; 
        while (std::getline(input, line)) {
          line = StringUtil::trim(line);
          if (line == "") {
            break;
          }
          peak_num++;
          line_list.push_back(line);
        }
        EnvelopePtr env_ptr(new Envelope(peak_num -1, line_list));
        env_ptrs_.push_back(env_ptr);        
      }
      input.close();
      initBaseMassIdx();
    }

    
// Initialize the mappping array between base peak masses and monoisotopic masses
void EnvDist::initBaseMassIdx() {
  base_mass_idxes_.resize(entry_num_, -1);
  for (int i = 0; i < entry_num_; i++) {
    double mz = env_ptrs_[i]->getReferMz();
    int idx = (int) (mz / mass_interval_);
    if (idx >= 0 && idx < entry_num_) {
      base_mass_idxes_[idx] = i;
    }
  }
  base_mass_idxes_[0] = 0;
  for (int i = 1; i < entry_num_; i++) {
    if (base_mass_idxes_[i] < 0) {
      base_mass_idxes_[i] = base_mass_idxes_[i - 1];
    }
  }
}

// Get a distribution envelope based on its monoisotopic mass.
EnvelopePtr EnvDist::getEnvByMonoMass(double mass) {
  int idx = (int) (mass / mass_interval_);
  if (idx < 0) {
    LOG_ERROR("Invalid mass");
    std::terminate();
  } else if (idx >= entry_num_) {
    LOG_DEBUG("mass out of bound");
    return nullptr;
  }
  return env_ptrs_[idx];
}

// Get a distribution envelope based on its mass of base peak.

EnvelopePtr EnvDist::getEnvByBaseMass(double mass) {
  int idx = (int) (mass / mass_interval_);
  if (idx < 0) {
    LOG_ERROR("Invalid mass");
    std::terminate();
  } else if (idx >= entry_num_) {
    LOG_DEBUG("Mass out of bound");
    return nullptr;
  }
  return env_ptrs_[base_mass_idxes_[idx]];
}

}
