#include "base/logger.hpp"
#include "feature/env_base.hpp" 

namespace prot {

EnvBase::EnvBase(std::string file_name, int entry_num, 
                 double mass_interval):
    entry_num_(entry_num),
    mass_interval_(mass_interval) {
      std::ifstream input;
      input.open(file_name.c_str(), std::ios::in);
      if (!input.is_open()) {
        LOG_ERROR( "env file  " << file_name << " does not exist.");
        throw "fasta file does not exist.";
      }

      for (int i = 0; i < entry_num_; i++) {
        int peak_num = 0;
        std::string line;
        std::string str;
        while (std::getline(input, line)) {
          if (line == "") {
            break;
          }
          peak_num++;
          str = str + line + "\n";
        }
        envs_[i] = EnvelopePtr(new Envelope(peak_num - 1, str));
      }
      input.close();
      initBaseMassIdx();
    }

void EnvBase::initBaseMassIdx() {
  base_mass_idxes_.resize(entry_num_, -1);
  for (int i = 0; i < entry_num_; i++) {
    double mz = envs_[i]->getReferMz();
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
EnvelopePtr EnvBase::getEnvByMonoMass(double mass) {
  int idx = (int) (mass / mass_interval_);
  if (idx < 0) {
    LOG_ERROR("Invalid mass");
    exit(1);
  } else if (idx >= entry_num_) {
    LOG_DEBUG("mass out of bound");
    return nullptr;
  }
  return envs_[idx];
}

// Get a distribution envelope based on its mass of base peak.
EnvelopePtr EnvBase::getEnvByBaseMass(double mass) {
  int idx = (int) (mass / mass_interval_);
  if (idx < 0) {
    LOG_ERROR("Invalid mass");
    exit(1);
  } else if (idx >= entry_num_) {
    LOG_DEBUG("Mass out of bound");
    return nullptr;
  }
  return envs_[base_mass_idxes_[idx]];
}

}
