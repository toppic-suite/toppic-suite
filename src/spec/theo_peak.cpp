#include <algorithm>
#include "spec/theo_peak.hpp"

namespace prot {

TheoPeak::TheoPeak(IonPtr ion_ptr,double unmod_mass,
                   double shift):
    Peak(unmod_mass + shift, 1.0),
    ion_ptr_(ion_ptr),
    shift_(shift) {
    }

} /* namespace prot */
