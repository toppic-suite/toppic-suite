#ifndef PROT_PRSM_PRSM_STR_HPP_
#define PROT_PRSM_PRSM_STR_HPP_

#include <memory>
#include <vector>
#include <string>

namespace prot {

class PrsmStr;
typedef std::shared_ptr<PrsmStr> PrsmStrPtr;

class PrsmStr {
 public:
  PrsmStr(const std::vector<std::string> &str_vec);

  std::vector<std::string> getStrVec() {return str_vec_;}

  int getSpectrumId() {return spectrum_id_;}

  int getDbSeqId() {return db_seq_id_;}

  std::string getDbSeqName() {return db_seq_name_;}
  
  double getMatchFragNum() {return match_frag_num_;}

  double getEValue() {return e_value_;}

  double getFdr() {return fdr_;}

  void setId(int id);

  void setFdr(double fdr);

  //comparison 
  static bool cmpEValueInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getEValue() < b->getEValue();
  }

  static bool cmpMatchFragmentDec(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getMatchFragNum() > b->getMatchFragNum();
  }

  static bool cmpSpectrumIdInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getSpectrumId() < b->getSpectrumId();
  }

 private:
  std::vector<std::string> str_vec_;

  int spectrum_id_;

  int db_seq_id_;

  std::string db_seq_name_;

  double match_frag_num_;

  double e_value_;

  double fdr_;
};

typedef std::vector<PrsmStrPtr> PrsmStrPtrVec;


}

#endif

