/*
 * comp_shift_hi_mem.hpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#ifndef COMP_SHIFT_HI_MEM_HPP_
#define COMP_SHIFT_HI_MEM_HPP_

#include "proteoform.hpp"
#include "ptm_fast_filter_mng.hpp"

namespace prot {

class CompShiftHiMem {
public:
	CompShiftHiMem(ProteoformPtrVec seqs,PtmFastFilterMngPtr mng);
	std::vector<std::vector<int>> compConvolution(std::vector<int> masses,int bgn_pos,int num);
	std::vector<std::vector<int>> compConvolution(std::vector<int> masses,std::vector<int> errors,int bgn_pos,int num);
	std::vector<std::vector<int>> getShiftScores(std::vector<short> scores,int num);
private:

	int shift_array_len_;
	int scale_;
	std::vector<int> seq_begins_;
	int seq_total_len_;
	std::vector<int> pos_seq_ids_;
	std::vector<int> index_begins_;
	std::vector<int> index_ends_;
	std::vector<int> indexes_;

	void initSeqBeginEnds(ProteoformPtrVec seqs);
	void initIndexex(ProteoformPtrVec seqs);
	void updateCnt(ProteoformPtr seq,std::vector<int> cnt);

};

} /* namespace prot */

#endif /* COMP_SHIFT_HI_MEM_HPP_ */
