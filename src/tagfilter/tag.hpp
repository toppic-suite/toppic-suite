
#ifndef PROT_TAG_HPP
#define PROT_TAG_HPP

#include <string>

#include "base/logger.hpp"

namespace prot {

class tag {

 public:
  tag (const std::string & X, const std::string & Y,
       double mass, bool ordered): X(X), Y(Y), 
    mass(mass), tolerance(200), ordered(ordered){}

  double getMinMass() {return mass - tolerance;}

  double getMaxMass() {return mass + tolerance;}

 private:
  std::string X, Y;
  double mass, tolerance;
  bool ordered;

};

}

#endif
