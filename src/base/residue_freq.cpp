#include "base/logger.hpp"
#include "base/residue_freq.hpp"

namespace prot {

ResidueFreq::ResidueFreq(AcidPtr acid_ptr, PtmPtr ptm_ptr, 
                         double freq): 
    Residue (acid_ptr, ptm_ptr),
    freq_(freq) {
    }

}
