#include "ion_type.hpp"

namespace prot {

IonType::IonType(std::string name, bool n_term, double shift) {
  name_ = name;
  n_term_ = n_term;
  shift_ = shift;
}


}
