#ifndef PROT_BASE_ALIGN_TYPE_HPP_
#define PROT_BASE_ALIGN_TYPE_HPP_

#include <string>
#include <memory>

namespace prot {

class AlignType;
typedef std::shared_ptr<AlignType> AlignTypePtr;

class AlignType {
 public: 
  static AlignTypePtr COMPLETE;
  static AlignTypePtr PREFIX;
  static AlignTypePtr SUFFIX;
  static AlignTypePtr INTERNAL;

  AlignType(const std::string &name, int id);

  std::string getName() {return name_;}

  int getId() {return id_;}

 private:
  std::string name_;
  int id_;
};


}

#endif
