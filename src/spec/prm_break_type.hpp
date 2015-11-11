#ifndef PROT_BASE_PRM_BREAK_TYPE_HPP_
#define PROT_BASE_PRM_BREAK_TYPE_HPP_

#include <memory>

namespace prot {

class PrmBreakType;
typedef std::shared_ptr<PrmBreakType> PrmBreakTypePtr;

class PrmBreakType {
 public:
  static const PrmBreakTypePtr NONE;
  static const PrmBreakTypePtr N_TERM;
  static const PrmBreakTypePtr C_TERM;
  static const PrmBreakTypePtr BOTH;

  std::string getName() {return name_;}

 private:
  std::string name_;
  PrmBreakType(std::string name): name_(name) {};
};

}

#endif

