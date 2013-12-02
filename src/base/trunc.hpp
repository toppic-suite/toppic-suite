#ifndef PROT_TRUNC_HPP_
#define PROT_TRUNC_HPP_

#include <string>
#include "base/acid.hpp"

namespace prot {

class Trunc {

 public:
  Trunc(std::string name, int trunc_len, 
           AcidPtrVec &acid_list, std::string str);
  std::string getName() {return name_;}
	int getTruncLen() {return trunc_len_;}
	double getShift() {return shift_;}

 private:
  std::string name_;
	int trunc_len_;
  AcidPtrVec acid_str_;
	double shift_;
};

typedef std::shared_ptr<Trunc> TruncPtr;
typedef std::vector<TruncPtr> TruncPtrVec;

TruncPtrVec getTruncPtrVecInstance(AcidPtrVec &acid_list, 
                                   const char* file_name);

TruncPtr getTruncPtrByName(TruncPtrVec &trunc_list, 
                           const std::string &name);
}

#endif
