#ifndef PROT_SEGMENT_HPP_
#define PROT_SEGMENT_HPP_

#include "base/proteoform.hpp"
#include "base/trunc.hpp"

namespace prot {

class Segment {
 public:
  Segment(int left_bp_pos, int right_bp_pos, double n_shift, double c_shift) {
    left_bp_pos_ = left_bp_pos;
    right_bp_pos_ = right_bp_pos;
    pep_n_term_shift_ = n_shift;
    pep_c_term_shift_ = c_shift;
  }

  int getLeftBpPos () {return left_bp_pos_;}

  int getRightBpPos() {return right_bp_pos_;}

  double getPepNTermShift() {return pep_n_term_shift_;}

  double getPepCTermShift() {return pep_c_term_shift_;}

 private:
  // segment begin and end are based on break_points
  int left_bp_pos_;
  int right_bp_pos_;
  double pep_n_term_shift_; 
  double pep_c_term_shift_;
};

typedef std::shared_ptr<Segment> SegmentPtr;
typedef std::vector<SegmentPtr> SegmentPtrVec;

TruncPtr findProtTermTrunc(TruncPtrVec truncs,int trunc_len,ResSeqPtr resseq){
	for(unsigned int i=0;i<truncs.size();i++){
		if(truncs[i]->isSameTrunc(trunc_len,resseq)){
			return truncs[i];
		}
	}
	return nullptr;
};

TruncPtr findProtNTermTrunc(ResSeqPtr seq,int trunc_len,TruncPtrVec allowed_trunc){
//	ResSeqPtr resseq = seq->getResSeqPtr();
	return findProtTermTrunc(allowed_trunc,trunc_len,seq);
};
TruncPtr findProtCTermTrunc(ResSeqPtr seq,int last_res_pos,TruncPtrVec allowed_trunc){
//	int trunc_len = seq->getResSeqPtr()->getLen()-1-last_res_pos;
//	ResSeqPtr resseq = seq->getResSeqPtr();
	int trunc_len = seq->getLen()-1-last_res_pos;
	return  findProtTermTrunc(allowed_trunc,trunc_len,seq);
};

bool isAlignPrefix(TruncPtr n_trunc,double pep_n_term_shift,double threshold){
	if(n_trunc != nullptr && pep_n_term_shift <= threshold){
		return true;
	}
	return false;
};

bool isAlignSuffix(TruncPtr c_trunc,double pep_c_term_shift,double threshold){
	if(c_trunc != nullptr && pep_c_term_shift <= threshold){
		return true;
	}
	return false;
}


}
#endif
