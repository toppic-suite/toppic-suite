#ifndef PROT_SIMPLE_PRSM_STR_HPP_
#define PROT_SIMPLE_PRSM_STR_HPP_

#include <memory>
#include <vector>
#include <string>

namespace prot {

class SimplePrsmStr {
 public:
  SimplePrsmStr(const std::vector<std::string> &str_vec);

  std::vector<std::string> getStrVec() {return str_vec_;}

  int getSpectrumId() {return spectrum_id_;}
  
  double getScore() {return score_;}

 private:
  std::vector<std::string> str_vec_;

  int spectrum_id_;

  double score_;

};

typedef std::shared_ptr<SimplePrsmStr> SimplePrsmStrPtr;
typedef std::vector<SimplePrsmStrPtr> SimplePrsmStrPtrVec;

inline bool simplePrsmStrScoreDown(const SimplePrsmStrPtr &a, const SimplePrsmStrPtr &b) {
  return a->getScore() > b->getScore();
}

std::string getValueStr(std::string line);

std::string getXmlLine(const std::vector<std::string> &str_vec,
                       const std::string &property);

}

#endif

