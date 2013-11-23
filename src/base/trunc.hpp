#ifndef PROT_TRUNC_H_
#define PROT_TRUNC_H_

#include <string>
#include "acid.hpp"

namespace prot {

class Trunc {

 public:
  Trunc(std::string name, int trunc_len, 
           AcidPtrVec &acid_ptr_vec, std::string acid_str);
  std::string getName() {return name_;}
	int getTruncLen() {return trunc_len_;}
	double getShift() {return shift_;}

 private:
  std::string name_;
	int trunc_len_;
  AcidPtrVec acid_ptr_str_;
	double shift_;
};

typedef std::shared_ptr<Trunc> TruncPtr;
typedef std::vector<TruncPtr> TruncPtrVec;

}

#endif
