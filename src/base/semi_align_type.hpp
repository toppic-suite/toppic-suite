#ifndef PROT_SEMI_ALIGN_TYPE_HPP_
#define PROT_SEMI_ALIGN_TYPE_HPP_

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

class SemiAlignTypeFactory {
 private:
  static SemiAlignTypePtr semi_align_type_complete_;
  static SemiAlignTypePtr semi_align_type_prefix_;
  static SemiAlignTypePtr semi_align_type_suffix_;
  static SemiAlignTypePtr semi_align_type_internal_;

 public:
  static SemiAlignTypePtr getCompletePtr() {return semi_align_type_complete_;}
  static SemiAlignTypePtr getPrefixPtr() {return semi_align_type_prefix_;}
  static SemiAlignTypePtr getSuffixPtr() {return semi_align_type_suffix_;}
  static SemiAlignTypePtr getInternalPtr() {return semi_align_type_internal_;}

  static std::string convertSemiAlignTypeToString (int i);
};

}

#endif
