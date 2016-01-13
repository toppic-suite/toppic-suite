#ifndef PROT_BASE_RM_BREAK_TYPE_HPP_
#define PROT_BASE_RM_BREAK_TYPE_HPP_

#include <memory>
#include <string>

namespace prot {

// residue mass break type
class RmBreakType;
typedef std::shared_ptr<RmBreakType> RmBreakTypePtr;

class RmBreakType {
 public:
  static const RmBreakTypePtr NONE;
  static const RmBreakTypePtr N_TERM;
  static const RmBreakTypePtr C_TERM;
  static const RmBreakTypePtr BOTH;

  std::string getName() {return name_;}

 private:
  std::string name_;
  RmBreakType(std::string name): name_(name) {};
};

}

#endif

