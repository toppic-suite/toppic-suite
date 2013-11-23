
#ifndef PROT_EXP_TRUNC_H_
#define PROT_EXP_TRUNC_H_

#include <string>
#include "acid.hpp"

namespace prot {

class ExpTrunc {

 public:
  ExpTrunc(int trunc_len, AcidPtrVec &acid_ptr_vec, std::string acid_str);
	int getTruncLen() {return trunc_len_;}
	double getShift() {return shift_;}

 private:
	int trunc_len_;
  AcidPtrVec acid_ptr_str_;
	double shift_;
};

typedef std::shared_ptr<ExpTrunc> ExpTruncPtr;
typedef std::vector<ExpTruncPtr> ExpTruncPtrVec;

}

#endif
