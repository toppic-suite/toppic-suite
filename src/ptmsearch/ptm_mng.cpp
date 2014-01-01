/*
 * ptm_mng.cpp
 *
 *  Created on: Dec 19, 2013
 *      Author: xunlikun
 */

#include "ptm_mng.hpp"

namespace prot {
PtmMng::PtmMng(){
	allow_prot_N_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"NONE"));
	allow_prot_N_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"ACETYLATION"));
	allow_prot_N_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"NME"));
	allow_prot_N_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"NME_ACETYLATION"));

	allow_prot_C_mods_.push_back(prot::getProtModPtrByName(base_data_->getProtModPtrVec(),"NONE"));

	allow_prot_N_truncs_.push_back(prot::getTruncPtrByName(base_data_->getTruncPtrVec(),"NONE"));
	allow_prot_N_truncs_.push_back(prot::getTruncPtrByName(base_data_->getTruncPtrVec(),"NME"));

	allow_prot_C_truncs_.push_back(prot::getTruncPtrByName(base_data_->getTruncPtrVec(),"NONE"));

	allow_pep_N_mods_.push_back(prot::findEmptyPtmPtr(base_data_->getPtmPtrVec()));

	allow_pep_C_mods_.push_back(prot::findEmptyPtmPtr(base_data_->getPtmPtrVec()));
}
} /* namespace prot */
