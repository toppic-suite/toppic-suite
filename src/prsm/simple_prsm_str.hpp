#ifndef PROT_PRSM_SIMPLE_PRSM_STR_HPP_
#define PROT_PRSM_SIMPLE_PRSM_STR_HPP_

#include <memory>
#include <vector>
#include <string>

namespace prot {

class SimplePrsmStr;

typedef std::shared_ptr<SimplePrsmStr> SimplePrsmStrPtr;

class SimplePrsmStr {
 public:
  SimplePrsmStr(const std::vector<std::string> &str_vec);

  std::vector<std::string> getStrVec() {return str_vec_;}

  int getSpectrumId() {return spectrum_id_;}
  
  double getScore() {return score_;}

  static bool cmpScoreDec(const SimplePrsmStrPtr &a, const SimplePrsmStrPtr &b) {
    return a->getScore() > b->getScore();
  }

 private:
  std::vector<std::string> str_vec_;

  int spectrum_id_;

  double score_;

};

typedef std::vector<SimplePrsmStrPtr> SimplePrsmStrPtrVec;

}

#endif

