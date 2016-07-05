#ifndef PROT_FEATURE_ENV_DIST_HPP_
#define PROT_FEATURE_ENV_DIST_HPP_

#include <memory>
#include <vector>

#include "feature/envelope.hpp"

namespace prot {

class EnvDist {
 public:
  EnvDist(std::string file_name, int entry_num, double mass_interval);

  EnvelopePtr getEnvByMonoMass(double mass);

  EnvelopePtr getEnvByBaseMass(double mass);

 private:
  int entry_num_;
  // the mass interval between two neighboring entries 
  double mass_interval_;
  // the list of distribution envelopes 
  EnvelopePtrVec env_ptrs_;
  // mapping distribution envelopes to the mass value of base peak 
  std::vector<int> base_mass_idxes_;

  void initBaseMassIdx();
};

}

#endif
