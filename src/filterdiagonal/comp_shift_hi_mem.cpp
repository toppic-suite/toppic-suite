/*
 * comp_shift_hi_mem.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#include <time.h>
#include <log4cxx/logger.h>

#include "comp_shift_hi_mem.hpp"

namespace prot {

log4cxx::LoggerPtr compShiftHiMem_logger(log4cxx::Logger::getLogger("CompShiftHiMem"));

CompShiftHiMem::CompShiftHiMem(ProteoformPtrVec seqs,PtmFastFilterMngPtr mng,IonTypePtrVec ion_type_ptr_vec){
	scale_ = mng->ptm_fast_filter_scale_;
	LOG4CXX_DEBUG(compShiftHiMem_logger, "Scale"+scale_);
	LOG4CXX_DEBUG(compShiftHiMem_logger, "Sequence number"+seqs.size());
	shift_array_len_ = 20000 * scale_ + 2;
	initSeqBeginEnds(seqs);
	initIndexes(seqs,ion_type_ptr_vec);
	LOG4CXX_DEBUG(compShiftHiMem_logger, "shift_array_len_ ="+shift_array_len_);
	LOG4CXX_DEBUG(compShiftHiMem_logger, "seq_total_len_"+seq_total_len_);
	LOG4CXX_DEBUG(compShiftHiMem_logger, "indexes.length"+indexes_.size());
}
CompShiftHiMem::~CompShiftHiMem(){
	delete pos_seq_ids_;
}
std::vector<std::vector<int>> CompShiftHiMem::compConvolution(std::vector<int> masses,int bgn_pos,int num){
	std::vector<short> scores;
	for(int i=0;i<seq_total_len_ ;i++){
		scores.push_back(0);
	}

	int begin_index =0;
	int end_index = 0;
	unsigned int m =0;

//	time_t time_s = time(NULL);
//	long time = (long)time_s;

	for(unsigned int i =bgn_pos+1;i<masses.size();i++){
		m = masses[i]-masses[bgn_pos];
		if(m>=shift_array_len_){
			break;
		}
		if(m>0){
			begin_index = index_begins_[m-1];
		}
		else{
			// if m<0?
			begin_index = index_begins_[m];
		}
		end_index = index_ends_[m+1];
		//logger
		for(int j = begin_index;j<end_index;j++){
			scores[indexes_[j]];
		}
	}
	//logger
	return getShiftScores(scores,num);
}
std::vector<std::vector<int>> CompShiftHiMem::compConvolution(std::vector<int> masses,std::vector<int> errors,int bgn_pos,int num){
	std::vector<short> scores;
	for(int i=0;i<seq_total_len_ ;i++){
		scores.push_back(0);
	}

	int begin_index =0;
	int end_index = 0;
	int m =0;
	int e =0;

//	time_t time_s = time(NULL);
//	long time = (long)time_s;

	for(unsigned int i =bgn_pos+1;i<masses.size();i++){
		m = masses[i]-masses[bgn_pos];
		int left = m-errors[i]-e;
		if(left < 0){
			left=0;
		}
		unsigned int right = m+errors[i]+e;
		if(right>=shift_array_len_-1){
			break;
		}
		begin_index = index_begins_[left];
		end_index= index_ends_[right];
		//logger
		for(int j=begin_index;j<=end_index;i++){
			scores[indexes_[j]]++;
		}
	}
	//logger
	return getShiftScores(scores,num);
}
std::vector<std::vector<int>> CompShiftHiMem::getShiftScores(std::vector<short> scores,int num){
//	time_t time_s = time(NULL);
//	long time = (long)time_s;

	//get top pos
	std::vector<short> top_scores;
	std::vector<int> top_position;
	for(int i =0; i<num;i++){
		top_scores.push_back(0);
		top_position.push_back(-1);
	}
	int last_scr = top_scores[num-1];
	for(int i=0;i<seq_total_len_;i++){
		if(scores[i]<last_scr){
			continue;
		}
		for(int j=num-2;j>=0;j--){
			if(scores[i]>top_scores[j]){
				top_scores[j+1] = top_scores[j];
				top_position[j+1] = top_position[j];
			}
			else{
				top_scores[j+1] = scores[i];
				top_position[j+1] = i;
				break;
			}
		}
		if(scores[i]>top_scores[0]){
			top_scores[0] = scores[i];
			top_position[0]=i;
		}
		last_scr = top_scores[num-1];
	}
	//logger
	int output_num =0;
	for(int i=0;i<num;i++){
		if(top_position[i]>=0){
			output_num++;
		}
	}
	std::vector<std::vector<int>> result;
	for(int i=0;i<output_num;i++){
		std::vector<int> temp;
		temp.push_back(pos_seq_ids_[top_position[i]]);
		temp.push_back(top_scores[i]);
		result.push_back(temp);
	}
	return result;
}

void CompShiftHiMem::initSeqBeginEnds(ProteoformPtrVec seqs){
	std::vector<int> seq_ends;
	int pnt = 0;
	for(unsigned int i=0;i<seqs.size();i++){
		seq_begins_.push_back(pnt);
		seq_ends.push_back(pnt+seqs[i]->getResSeqPtr()->getLen()-1);
		pnt+= seqs[i]->getResSeqPtr()->getLen();
	}
	seq_total_len_ = pnt;
	pos_seq_ids_ = new int[seq_total_len_];
	for(unsigned int i =0;i<seqs.size();i++){
		for(int j= seq_begins_[i];j<=seq_ends[i];j++){
			pos_seq_ids_[j] = i;
		}
	}
}
void CompShiftHiMem::initIndexes(ProteoformPtrVec seqs,IonTypePtrVec ion_type_ptr_vec){

	std::vector<int> cnt;
	std::vector<int> index_pnt;
	for(unsigned int i =0;i< shift_array_len_;i++){
		cnt.push_back(0);
		index_begins_.push_back(0);
		index_ends_.push_back(0);
		index_pnt.push_back(0);
	}

	for(unsigned int i =0;i<seqs.size();i++){
		CompShiftHiMem::updateCnt(seqs[i],cnt,ion_type_ptr_vec);
	}
	int pnt = 0;
	for(unsigned int i=0;i<cnt.size();i++){
		index_begins_[i] = pnt;
		index_pnt[i] = pnt;
		index_ends_[i] = pnt;
		pnt+=cnt[i];
	}

	for(int i=0;i<pnt;i++){
		indexes_.push_back(0);
	}

	for(unsigned int i=0;i<seqs.size();i++){
		if(i/10000*1000 == i){
			LOG4CXX_DEBUG(compShiftHiMem_logger, "preprocessing seq "+i);
		}
		std::vector<int> mass = seqs[i]->getBpSpecPtr()->getScaledMass(scale_,prot::getIonTypePtrByName(ion_type_ptr_vec,"B"));
		unsigned int bgn = 0;
		unsigned int end =0;
		unsigned int diff =0;
		while(bgn <mass.size()){
			while(end < mass.size()){
				diff = mass[end] - mass[bgn];
				if(diff >= shift_array_len_){
					break;
				}
				else{
					indexes_[index_pnt[diff]] =seq_begins_[i];
					index_pnt[diff]++;
				}
				end++;
			}
			bgn++;
		}
	}
}
void CompShiftHiMem::updateCnt(ProteoformPtr seq,std::vector<int> cnt,IonTypePtrVec ion_type_ptr_vec){

	std::vector<int> mass = seq->getBpSpecPtr()->getScaledMass(scale_,prot::getIonTypePtrByName(ion_type_ptr_vec,"B"));
	unsigned int bgn =0;
	unsigned int end =0;
	unsigned int diff=0;
	while(bgn < mass.size()){
		end = bgn + 1;
		while (end < mass.size()){
			diff = mass[end]-mass[bgn];
			if(diff >= cnt.size()){
				break;
			}
			else{
				cnt[diff]++;
			}
			end ++;
		}
		bgn++;
	}
}


} /* namespace prot */
