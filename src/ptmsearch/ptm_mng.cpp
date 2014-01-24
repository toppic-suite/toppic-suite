/*
 * ptm_mng.cpp
 *
 *  Created on: Dec 19, 2013
 *      Author: xunlikun
 */

#include "ptm_mng.hpp"

namespace prot {
PtmMng::PtmMng(std::string config_file_name){
	base_data_=BaseDataPtr(new BaseData(config_file_name));
	allow_prot_N_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"NONE"));
	allow_prot_N_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"ACETYLATION"));
	allow_prot_N_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"NME"));
	allow_prot_N_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"NME_ACETYLATION"));

	allow_prot_C_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"NONE"));

	allow_prot_N_truncs_.push_back(TruncFactory::getBaseTruncPtrByName("NONE"));
	allow_prot_N_truncs_.push_back(TruncFactory::getBaseTruncPtrByName("NME"));

	allow_prot_C_truncs_.push_back(TruncFactory::getBaseTruncPtrByName("NONE"));

	allow_pep_N_mods_.push_back(PtmFactory::findEmptyPtmPtr());

	allow_pep_C_mods_.push_back(PtmFactory::findEmptyPtmPtr());
}
} /* namespace prot */
