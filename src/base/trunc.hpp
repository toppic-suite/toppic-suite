#ifndef PROT_TRUNC_HPP_
#define PROT_TRUNC_HPP_

#include <string>
#include "base/acid.hpp"
#include "base/residue_seq.hpp"

namespace prot {

class Trunc {

 public:
  Trunc(std::string name, int trunc_len, std::string str);

  std::string getName() {return name_;}

  int getTruncLen() {return trunc_len_;}

  AcidPtrVec getAcidPtrVec() {return acid_str_;}

  double getShift() {return shift_;}

  void appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  bool isSameTrunc(int len, ResSeqPtr resseq);

 private:
  std::string name_;
	int trunc_len_;
  AcidPtrVec acid_str_;
	double shift_;
};

typedef std::shared_ptr<Trunc> TruncPtr;
typedef std::vector<TruncPtr> TruncPtrVec;

TruncPtr findProtNTermTrunc(ResSeqPtr seq,int trunc_len,TruncPtrVec allowed_trunc);

TruncPtr findProtCTermTrunc(ResSeqPtr seq,int last_res_pos,TruncPtrVec allowed_trunc);

bool isAlignPrefix(TruncPtr n_trunc,double pep_n_term_shift,double threshold);

bool isAlignSuffix(TruncPtr c_trunc,double pep_c_term_shift,double threshold);

/* trunc factory */
class TruncFactory {
 public:
  static void initFactory(const std::string &file_name);
  static TruncPtrVec& getBaseTruncPtrVec() {return trunc_ptr_vec_;}
  static TruncPtr getBaseTruncPtrByName(const std::string &name);

 private:
  static TruncPtrVec trunc_ptr_vec_;
};

}

#endif
