/*
 * ps_align.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: xunlikun
 */

#include <ptmsearch/ps_align.hpp>
#include <float.h>
#include "prsm/sem_align_type.hpp"

namespace prot {
PSAlign::PSAlign(){};
PSAlign::PSAlign(std::vector<double> sp_masses,std::vector<double> seq_masses,BasicDiagPairDiagPtrVec diagonals,PtmMngPtr mng){
	mng_ =mng;
	sp_masses_=sp_masses;
	seq_masses_=seq_masses;
	diagonals_ = diagonals;
}
void PSAlign::compute(int align_type){
	initDPPair();
	dp(align_type);
	backtrace(align_type);
}
void PSAlign::initDPPair(){
	for(int i=0;i<diagonals_.size();i++){
		segment_bgn_pairs_.push_back(nullptr);
		segment_end_pairs_.push_back(nullptr);
		DPPairPtrVec temp_dppair;
		dp_2d_pairs_.push_back(temp_dppair);
		for(int j =0;j<diagonals_[i]->size();j++){
			int x = diagonals_[i]->getDiagPair(j)->getX();
			int y = diagonals_[i]->getDiagPair(j)->getY();
			double score = diagonals_[i]->getDiagPair(j)->getScore();
			double diff = diagonals_[i]->getDiagPair(j)->getDiff();
			dp_2d_pairs_[i].push_back(DPPairPtr(new DPPair(x,y,score,diff,j,mng_->n_unknown_shift_,diagonals_[i]->getHeader())));
		}
		segment_bgn_pairs_[i] = dp_2d_pairs_[i][0];
		segment_end_pairs_ = dp_2d_pairs_[i][diagonals_[i]->size()-1];
	}

	first_pair_ =DPPairPtr( new DPPair(-1,-1,0,0,-1,mng_->n_unknown_shift_,nullptr));
	first_pair_->setDiagPrev(nullptr);
	dp_pairs_.push_back(first_pair_);
	for(int i=0;i<dp_2d_pairs_.size();i++){
		for(int j=0;j<dp_2d_pairs_[i].size();j++){
			dp_pairs_.push_back(dp_2d_pairs_[i][j]);
			if(j>0){
				dp_2d_pairs_[i][j]->setDiagPrev(dp_2d_pairs_[i][j-1]);
			}
		}
	}

	double diff = sp_masses_[sp_masses_.size()-1]-seq_masses_[seq_masses_.size()-1];
	last_pair_ = DPPairPtr( new DPPair(sp_masses_.size(),seq_masses_.size(),0,diff,-1,mng_->n_unknown_shift_,nullptr));
	last_pair_ ->setDiagPrev(nullptr);
	dp_pairs_.push_back(last_pair_);

}
void PSAlign::dpPrep(){
	std::sort(dp_pairs_.begin(),dp_pairs_.end(),prot::comparePairUp);
	for(int s=0;s<mng_->n_unknown_shift_+1;s++){
		dp_pairs_[0]->updateTable(s,0,PATH_TYPE_NULL,nullptr);
	}
}

DPPairPtr PSAlign::getTruncPre(DPPairPtr cur_pair,int s,int type){
	DPPairPtr trunc_prev;
	if(cur_pair == last_pair_){
		double trunc_score = -DBL_MAX;
		for(int i=0;i<segment_end_pairs_.size();i++){
			DPPairPtr prev_pair = segment_end_pairs_[i];
			if(type == SEMI_ALIGN_TYPE_COMPLETE || type == SEMI_ALIGN_TYPE_SUFFIX){
				if(prev_pair->getDiagonalHeader()->isAllowProtCMod()&& prev_pair->getSrc(s)>trunc_score){
					trunc_prev = prev_pair;
					trunc_score = prev_pair->getSrc(s);
				}
			}
			else{
				if(prev_pair->getDiagonalHeader()->isAllowPepCMod()&& prev_pair->getSrc(s)>trunc_score){
					trunc_prev = prev_pair;
					trunc_score = prev_pair->getSrc(s);
				}
			}
		}
	}
	else{
		if(cur_pair->getDiagOrder() == 0){
			if(type == SEMI_ALIGN_TYPE_COMPLETE || type == SEMI_ALIGN_TYPE_PREFIX){
				if(cur_pair->getDiagonalHeader()->isAllowProtNMod()){
					trunc_prev = first_pair_;
				}
			}
			else{
				if(cur_pair->getDiagonalHeader()->isAllowPepNMod()){
					trunc_prev = first_pair_;
				}
			}
		}
	}
	return trunc_prev;
}
DPPairPtr PSAlign::getShiftPre(DPPairPtr cur_pair,int p,int s,int type){
	int cur_x = cur_pair->getX();
	int cur_y = cur_pair->getY();
	DPPairPtr shift_prev = nullptr;
	double shift_score = -DBL_MAX;
	if(s>=1){
		if(cur_pair == last_pair_){
			for(int q=0;q<p;q++){
				DPPairPtr prev_pair = dp_pairs_[q];
				int prev_x = prev_pair->getX();
				int prev_y = prev_pair->getY();
				if(prev_x >= cur_x || prev_y >=cur_y ||prev_pair->getDiagonalHeader() == cur_pair->getDiagonalHeader()){
					continue;
				}
				if(type == SEMI_ALIGN_TYPE_COMPLETE || type == SEMI_ALIGN_TYPE_SUFFIX){
					continue;
				}
				if(prev_pair->getSrc(s-1)>shift_score){
					shift_prev = prev_pair;
					shift_score = dp_pairs_[q]->getSrc(s-1);
				}
			}
		}
		else{
			if(type == SEMI_ALIGN_TYPE_COMPLETE || type == SEMI_ALIGN_TYPE_PREFIX){
				if(first_pair_->getSrc(s-1)>shift_score&& cur_pair->getDiagonalHeader()->isAlignPrefix()){
					shift_prev = first_pair_;
					shift_score = first_pair_->getSrc(s-1);
				}
			}
			else{
				if(first_pair_->getSrc(s-1)>shift_score){
					shift_prev = first_pair_;
					shift_score = first_pair_->getSrc(s-1);
				}
			}
			for(int q=1;q<p;q++){
				DPPairPtr prev_pair = dp_pairs_[q];
				int prev_x = prev_pair->getX();
				int prev_y = prev_pair->getY();
				double gap = std::abs(prev_pair->getDiagonalHeader()->getProtNTermShift()-cur_pair->getDiagonalHeader()->getProtNTermShift());
				if(prev_x>=cur_x||prev_y>=cur_y||prev_pair->getDiagonalHeader()== cur_pair->getDiagonalHeader()||gap<=mng_->align_min_gap){
					continue;
				}
				double prev_score = prev_pair->getSrc(s-1);
				if(gap>mng_->large_shift_thresh){
					prev_score = prev_score - mng_->large_shift_panelty;
				}
				if(prev_score > shift_score){
					shift_prev = prev_pair;
					shift_score = prev_score;
				}
			}
		}
	}
	return shift_prev;
}
void PSAlign::dp(int align_type){
	dpPrep();
	for(int p =1;p<dp_pairs_.size();p++){
		for(int s =0;s<mng_->n_unknown_shift_;s++){
			DPPairPtr trunc_prev = getTruncPre(dp_pairs_[p],s,align_type);
			double trunc_score;
			if(trunc_prev == nullptr){
				trunc_score = -DBL_MAX;
			}
			else{
				trunc_score =trunc_prev->getSrc(s);
			}
			DPPairPtr diag_prev = dp_pairs_[p]->getDiagPrev();
			double diag_score;
			if(diag_prev!=nullptr){
				diag_score = diag_prev->getSrc(s);
			}
			else{
				diag_score = -DBL_MAX;
			}
			DPPairPtr shift_prev = getShiftPre(dp_pairs_[p],p,s,align_type);
			double shift_score;
			if(shift_prev == nullptr){
				shift_score= -DBL_MAX;
			}
			else{
				shift_score = shift_prev->getSrc(s-1);
			}

			double new_score = dp_pairs_[p]->getPairScore();
			if(trunc_score >= diag_score && trunc_score >= shift_score ){
				if(trunc_score == -DBL_MAX){
					dp_pairs_[p]->updateTable(s,-DBL_MAX,PATH_TYPE_NULL,nullptr);
				}
				else{
					dp_pairs_[p]->updateTable(s,trunc_score+new_score,PATH_TYPE_TRUNC,trunc_prev);
				}
			}
			else if(diag_score >= shift_score){
				dp_pairs_[p]->updateTable(s,diag_score+new_score,PATH_TYPE_DIAGONAL,diag_prev);
			}
			else{
				dp_pairs_[p]->updateTable(s,shift_score+new_score,PATH_TYPE_SHIFT,shift_prev);
			}
		}
	}
}
void PSAlign::backtrace(int align_type){
	for(int s=0;s<=mng_->n_unknown_shift_;s++){
		align_scores_.push_back(0.0);
		DiagonalHeaderPtrVec temp;
		backtrack_diagonals_.push_back(temp);
		backtrack_diagonals_[s] = backtrace(s,align_type);
	}

}
DiagonalHeaderPtrVec PSAlign::backtrace(int s,int type){
	DiagonalHeaderPtrVec list;
	DiagonalHeaderPtr cur_header;
	int cur_end = -1;
	int cur_bgn = -1;
	DPPairPtr p = last_pair_;

	align_scores_[s] = p->getSrc(s);

	if(p->getPre(s)==nullptr || p->getPre(s)==first_pair_){
		return list;
	}

	while(p!=nullptr){
		DPPairPtr pre = p->getPre(s);
		if(p==last_pair_){
			cur_header = pre->getDiagonalHeader();
			cur_end = pre->getY()-1;
		}
		else if(pre == first_pair_){
			cur_bgn = p->getY();
			list.push_back(prot::getShift(cur_header,cur_bgn,cur_end));
		}
		else{
			if(p->getType(s)==PATH_TYPE_SHIFT){
				cur_bgn=p->getY();
				list.push_back(prot::getShift(cur_header,cur_bgn,cur_end));
				cur_header = pre->getDiagonalHeader();
				cur_end = pre->getY()-1;
			}
		}
		if(p->getType(s)==PATH_TYPE_SHIFT){
			s--;
		}
		p=pre;
	}
	DiagonalHeaderPtrVec sub_seg;
	for(int i=0;i<list.size();i++){
		sub_seg.push_back(list[list.size()-1-i]);
	}
	return sub_seg;
}
} /* namespace prot */
