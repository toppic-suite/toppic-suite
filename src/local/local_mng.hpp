#ifndef PROT_LOCAL_MNG_HPP_
#define PROT_LOCAL_MNG_HPP_

#include "prsm/prsm_para.hpp"
#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"

namespace prot {

const double SCALE_FACTOR = 274.335215;
const int LEFT_SUP_LIMIT = 10;
const int RIGHT_SUP_LIMIT = 10;

class LocalMng {
  public:

    LocalMng(PrsmParaPtr prsm_para_ptr, const std::string& protect, const std::string& t,
             double max_ptm_mass, const std::string &input_file_ext,
             const std::string &output_file_ext);

    std::string input_file_ext_;
    std::string output_file_ext_;
    double thread_, theta_, beta_, weight_;
    PrsmParaPtr prsm_para_ptr_;
    bool cysteine_protected_;
    double min_mass_, max_ptm_mass_;
    double p1_, p2_;
};

typedef std::shared_ptr<LocalMng> LocalMngPtr;

inline LocalMng::LocalMng(PrsmParaPtr prsm_para_ptr, const std::string& protect,
                          const std::string& t, double max_ptm_mass,
                          const std::string& input_file_ext,
                          const std::string& output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;
    if (protect == "C0") {
        cysteine_protected_ = false;
    } else {
        cysteine_protected_ = true;
    }
    min_mass_ = prsm_para_ptr->getSpParaPtr()->getMinMass();
    thread_ = std::stod(t);
    max_ptm_mass_ = max_ptm_mass;
    weight_ = theta_ = 0.994;
    beta_ = 0.28;
    p1_ = p2_ = 0.0;
}

}

#endif /* PROT_LOCAL_MNG_HPP_ */
