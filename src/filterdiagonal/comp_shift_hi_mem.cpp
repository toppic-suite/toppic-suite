/*
 * comp_shift_hi_mem.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#include <log4cxx/logger.h>

#include "comp_shift_hi_mem.hpp"

namespace prot {

log4cxx::LoggerPtr compShiftHiMem_logger(log4cxx::Logger::getLogger("CompShiftHiMem"));

CompShiftHiMem::CompShiftHiMem(ProteoformPtrVec seqs,PtmFastFilterMngPtr mng){
	scale_ = mng->ptm_fast_filter_scale_;
	LOG4CXX_DEBUG(compShiftHiMem_logger, "Scale"+scale_);
	LOG4CXX_DEBUG(compShiftHiMem_logger, "Sequence number"+seqs.size());
	shift_array_len_ = 20000 * scale_ + 2;
// todo: implement init function
//	initSeqBeginEnds(seqs);
//	initIndexex(seqs);
	LOG4CXX_DEBUG(compShiftHiMem_logger, "shift_array_len_ ="+shift_array_len_);
	LOG4CXX_DEBUG(compShiftHiMem_logger, "seq_total_len_"+seq_total_len_);
	LOG4CXX_DEBUG(compShiftHiMem_logger, "indexes.length"+indexes_.size());
}
std::vector<std::vector<int>> CompShiftHiMem::compConvolution(std::vector<int> masses,int bgn_pos,int num){

}
std::vector<std::vector<int>> CompShiftHiMem::compConvolution(std::vector<int> masses,std::vector<int> errors,int bgn_pos,int num){

}
std::vector<std::vector<int>> CompShiftHiMem::getShiftScores(std::vector<short> scores,int num){

}

void initSeqBeginEnds(ProteoformPtrVec seqs){

}
void initIndexex(ProteoformPtrVec seqs){

}
void updateCnt(ProteoformPtr seq,std::vector<int> cnt){

}


} /* namespace prot */
