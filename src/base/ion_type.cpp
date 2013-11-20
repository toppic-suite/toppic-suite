
#include <map>
#include "ion_type.hpp"

namespace prot {


const IonTypePtr IonType::B = IonTypePtr(new IonType("B", true, 0));
const IonTypePtr IonType::Y = IonTypePtr(new IonType("Y", false, 18.0106));
const IonTypePtr IonType::C = IonTypePtr(new IonType("C", true, 17.0265));
const IonTypePtr IonType::Z_DOT = IonTypePtr(new IonType("Z_DOT", false, 1.9919));

static std::map<std::string, IonTypePtr> ion_type_ptr_map_ = {
  {"B", IonType::B}, {"Y",IonType::Y}, 
  {"C", IonType::C}, {"Z_DOT", IonType::Z_DOT}};

IonType::IonType(std::string name, bool n_term, double shift) {
  name_ = name;
  n_term_ = n_term;
  shift_ = shift;
}

/* Gets IonType by name */
IonTypePtr IonType::getIonTypePtrByName(std::string name) {
  return ion_type_ptr_map_[name];
}

}
