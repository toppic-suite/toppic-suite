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

  std::string getSeqName() {return seq_name_;}

  int getSpeciesId() {return species_id_;}

  double getMatchFragNum() {return match_frag_num_;}

  double getNormMatchFragNum() {return norm_match_frag_num_;}

  double getEValue() {return e_value_;}

  double getFdr() {return fdr_;}

  double getProteoformFdr() {return proteoform_fdr_;}

  void setId(int id);

  void setFdr(double fdr);

  void setProteoformFdr(double proteoform_fdr);

  //comparison 
  static bool cmpEValueInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getEValue() < b->getEValue();
  }

  static bool cmpMatchFragmentDec(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getMatchFragNum() > b->getMatchFragNum();
  }

  static bool cmpNormMatchFragmentDec(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getNormMatchFragNum() > b->getNormMatchFragNum();
  }

  static bool cmpSpectrumIdInc(const PrsmStrPtr &a, const PrsmStrPtr &b) {
    return a->getSpectrumId() < b->getSpectrumId();
  }

 private:
  std::vector<std::string> str_vec_;

  int spectrum_id_;

  std::string seq_name_;

  int species_id_;

  double match_frag_num_;

  /* a tempory variable for testing mass graph alignment */
  double norm_match_frag_num_;

  double e_value_;

  double fdr_;

  double proteoform_fdr_;
};

typedef std::vector<PrsmStrPtr> PrsmStrPtrVec;
typedef std::vector<PrsmStrPtrVec> PrsmStrPtr2D;

}

#endif

