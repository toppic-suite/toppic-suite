
#ifndef SUPPORT_PEAK_TYPE_HPP_
#define SUPPORT_PEAK_TYPE_HPP_

#include <memory>
#include <string>
#include <vector>

#include <xercesc/dom/DOM.hpp>

namespace prot {

#define SP_TYPE_N_TERM "N_TERM"

class SupportPeakType {
 public:
  SupportPeakType(int id, const std::string &name){
    id_=id;
    name_=name;
  }
  int getId(){return id_;}
  std::string getName(){return name_;}
 private:
  int id_;
  std::string name_;
};

typedef std::shared_ptr<SupportPeakType> SPTypePtr;
typedef std::vector<SPTypePtr> SPTypePtrVec;

/* support peak type factory */
class SPTypeFactory {
 public:
  static void initFactory(const std::string &file_name);
  static const SPTypePtrVec& getBaseSPTypePtrVec() {
    return sp_type_ptr_vec_;}

  static SPTypePtr getBaseSPTypePtrByName(const std::string &name);
  static SPTypePtr getBaseSPTypePtrById(int id);

  static SPTypePtr getSPTypePtr_N_TERM() {
    return getBaseSPTypePtrByName(SP_TYPE_N_TERM);
  }
 private:
  static SPTypePtrVec sp_type_ptr_vec_;
};

} /* namespace prot */

#endif /* SUPPORT_PEAK_TYPE_HPP_ */
