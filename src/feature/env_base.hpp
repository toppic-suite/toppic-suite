#ifndef PROT_FEATURE_ENVELOPE_BASE_HPP_
#define PROT_FEATURE_ENVELOPE_BASE_HPP_

#include <vector>

#include "feature/envelope.hpp"

namespace prot {

class EnvBase {
 public:
  EnvBase(std::string file_name_, int entry_num_, double mass_interval_);

  void initBaseMassIdx();

  EnvelopePtr getEnvByMonoMass(double mass);

  EnvelopePtr getEnvByBaseMass(double mass);

 private:
  // number of distribution entries 
  int entry_num_;
  // the mass interval between two neighboring entries 
  double mass_interval_;
  // the list of distribution envelopes 
  EnvelopePtrVec envs_;
  // mapping distribution envelopes to the mass value of base peak 
  std::vector<int> base_mass_idxes_;
};

typedef std::shared_ptr<EnvBase> EnvBasePtr;

typedef std::vector<EnvBasePtr> EnvBasePtrVec;

}

#endif
