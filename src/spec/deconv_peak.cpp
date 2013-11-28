
#include "deconv_peak.hpp"

namespace prot {

DeconvPeak::DeconvPeak (int id, double mono_mass, 
                        double intensity, int charge) 
    : Peak (mono_mass, intensity) {
      id_ = id;
      charge_ = charge;
      score_ = 1.0;
    }
}

