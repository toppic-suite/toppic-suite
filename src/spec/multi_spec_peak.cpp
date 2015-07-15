#include "spec/multi_spec_peak.hpp"

namespace prot {

MultiSpecPeak::MultiSpecPeak (double mono_mass, double intensity, int spec_id, 
                              int in_spec_peak_id, int base_type, double score)
    : Peak (mono_mass, intensity) {
      spec_id_ = spec_id;
      in_spec_peak_id_ = in_spec_peak_id;
      base_type_ = base_type;
      score_ = score;
    }
}

