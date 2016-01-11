#ifndef PROT_SPEC_BASE_PEAK_TYPE_HPP_
#define PROT_SPEC_BASE_PEAK_TYPE_HPP_

#include <memory>
#include <string>

namespace prot {

class BasePeakType;
typedef std::shared_ptr<BasePeakType> BasePeakTypePtr;

class BasePeakType {
 public:
  static const BasePeakTypePtr ORIGINAL;
  static const BasePeakTypePtr REVERSED;

  BasePeakType(std::string name) {name_=name;}

  std::string getName() {return name_;}

 private:
  std::string name_;
};


}

#endif

