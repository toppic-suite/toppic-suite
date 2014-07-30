#ifndef PROT_TRUNC_HPP_
#define PROT_TRUNC_HPP_

#include <string>
#include "base/acid.hpp"
#include "base/residue_seq.hpp"

namespace prot {

class Trunc {

 public:
  Trunc(const std::string &name, int trunc_len, 
        const std::string &str);

  const std::string& getName() {return name_;}

  int getTruncLen() {return trunc_len_;}

  const AcidPtrVec& getAcidPtrVec() {return acid_str_;}

  double getShift() {return shift_;}

  void appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  bool isSameTrunc(int len, ResSeqPtr res_seq_ptr);

  bool isValidTrunc(ResSeqPtr res_seq_ptr);

 private:
  std::string name_;
  int trunc_len_;
  AcidPtrVec acid_str_;
  double shift_;

};

typedef std::shared_ptr<Trunc> TruncPtr;
typedef std::vector<TruncPtr> TruncPtrVec;

/* trunc factory */
class TruncFactory {
 public:
  static void initFactory(const std::string &file_name);
  static const TruncPtrVec& getBaseTruncPtrVec() {return trunc_ptr_vec_;}
  static TruncPtr getBaseTruncPtrByName(const std::string &name);

 private:
  static TruncPtrVec trunc_ptr_vec_;
};

}

#endif
