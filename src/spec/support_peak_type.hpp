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

typedef std::shared_ptr<SupportPeakType> SupportPeakTypePtr;
typedef std::vector<SupportPeakTypePtr> SupportPeakTypePtrVec;

SupportPeakTypePtrVec getSupportPeakTypePtrVecInstance(const char* file_name);

SupportPeakTypePtr getSupportPeakTypePtrByName(SupportPeakTypePtrVec &support_peak_type_list,
                         const std::string &name);

SupportPeakTypePtr getSupportPeakTypePtrById(SupportPeakTypePtrVec &support_peak_type_list,
                         const int id);

} /* namespace prot */

#endif /* SUPPORT_PEAK_TYPE_HPP_ */
