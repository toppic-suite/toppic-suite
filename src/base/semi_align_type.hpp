#ifndef PROT_BASE_SEMI_ALIGN_TYPE_HPP_
#define PROT_BASE_SEMI_ALIGN_TYPE_HPP_

#include <string>
#include <memory>

namespace prot {

class SemiAlignType {
 public: 
  SemiAlignType(const std::string &name, int id);

  std::string getName() {return name_;}

  int getId() {return id_;}

 private:
  std::string name_;
  int id_;
};

typedef std::shared_ptr<SemiAlignType> SemiAlignTypePtr;

}

#endif
