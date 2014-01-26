/*
 * support_peak_type.hpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#ifndef SUPPORT_PEAK_TYPE_HPP_
#define SUPPORT_PEAK_TYPE_HPP_

#include <memory>
#include <string>
#include <vector>

#include <xercesc/dom/DOM.hpp>

namespace prot {

class SupportPeakType {
 public:
  SupportPeakType(int id,std::string name){
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
  static SPTypePtrVec& getBaseSPTypePtrVec() {
    return sp_type_ptr_vec_;}

  static SPTypePtr getBaseSPTypePtrByName(const std::string &name);
  static SPTypePtr getBaseSPTypePtrById(const int id);

 private:
  static SPTypePtrVec sp_type_ptr_vec_;
};

} /* namespace prot */

#endif /* SUPPORT_PEAK_TYPE_HPP_ */
