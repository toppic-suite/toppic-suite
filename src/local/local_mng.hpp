#ifndef PROT_LOCAL_MNG_HPP_
#define PROT_LOCAL_MNG_HPP_

#include "prsm/prsm_para.hpp"
#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"

namespace prot {

const int LEFT_SUP_LIMIT = 10;
const int RIGHT_SUP_LIMIT = 10;
const int DESC_MATCH_LIMIT = 5;

typedef std::pair<PtmPtr, PtmPtr> PtmPair;
typedef std::vector<PtmPair> PtmPairVec;

class LocalMng {
 public:

  LocalMng(PrsmParaPtr prsm_para_ptr, const std::string& local_threshold, 
           const std::string& residueModFileName, double max_ptm_mass,
           const std::string &input_file_ext, const std::string &output_file_ext);

  PrsmParaPtr prsm_para_ptr_;

  std::string residueModFileName_;
  std::string input_file_ext_;
  std::string output_file_ext_;

  double thread_, theta_, beta_;
  double min_mass_, max_ptm_mass_;
  double p1_, p2_;
};

typedef std::shared_ptr<LocalMng> LocalMngPtr;

inline LocalMng::LocalMng(PrsmParaPtr prsm_para_ptr, const std::string& local_threshold, 
                          const std::string& residueModFileName, double max_ptm_mass,
                          const std::string &input_file_ext, const std::string &output_file_ext) {

  prsm_para_ptr_ = prsm_para_ptr;
  input_file_ext_ = input_file_ext;
  output_file_ext_ = output_file_ext;
  residueModFileName_ = residueModFileName;
  min_mass_ = prsm_para_ptr->getSpParaPtr()->getMinMass();
  thread_ = std::stod(local_threshold);
  max_ptm_mass_ = max_ptm_mass;
  theta_ = 0.994;
  beta_ = 0.8;
  p1_ = 0.915258;
  p2_ = 21.1822;
}

}

#endif /* PROT_LOCAL_MNG_HPP_ */
