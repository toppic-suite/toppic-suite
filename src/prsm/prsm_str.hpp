#ifndef PROT_PRSM_STR_HPP_
#define PROT_PRSM_STR_HPP_

#include <memory>
#include <vector>
#include <string>

namespace prot {

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

 private:
  std::vector<std::string> str_vec_;

  int spectrum_id_;

  int db_seq_id_;

  std::string db_seq_name_;

  double match_frag_num_;

  double e_value_;

  double fdr_;
};

typedef std::shared_ptr<PrsmStr> PrsmStrPtr;
typedef std::vector<PrsmStrPtr> PrsmStrPtrVec;

inline bool prsmStrEValueUp(const PrsmStrPtr &a, const PrsmStrPtr &b) {
  return a->getEValue() < b->getEValue();
}

inline bool prsmStrMatchFragmentDown(const PrsmStrPtr &a, const PrsmStrPtr &b) {
  return a->getMatchFragNum() > b->getMatchFragNum();
}

inline bool prsmStrSpectrumIdUp(const PrsmStrPtr &a, const PrsmStrPtr &b) {
  return a->getSpectrumId() < b->getSpectrumId();
}

}

#endif

