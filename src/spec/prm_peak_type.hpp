/*
 * prm_peak_type.hpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#ifndef PROT_PRM_PEAK_TYPE_HPP_
#define PROT_PRM_PEAK_TYPE_HPP_

#include <memory>
#include <vector>

#include <xercesc/dom/DOM.hpp>

namespace prot {

class PrmPeakType {
public:
	PrmPeakType(int id,std::string name){
		id_=id;
		name_=name;
	}
	int getId(){return id_;}
	std::string getName(){return name_;}
private:
	int id_;
	std::string name_;
};

typedef std::shared_ptr<PrmPeakType> PrmPeakTypePtr;
typedef std::vector<PrmPeakTypePtr> PrmPeakTypePtrVec;

PrmPeakTypePtrVec getPrmPeakTypePtrVecInstance(const char* file_name);

PrmPeakTypePtr getPrmPeakTypePtrByName(PrmPeakTypePtrVec &prm_peak_type_list,
                         const std::string &name);
PrmPeakTypePtr getPrmPeakTypePtrById(PrmPeakTypePtrVec &prm_peak_type_list,
                         const int id );
} /* namespace prot */

#endif /* PRM_PEAK_TYPE_HPP_ */
