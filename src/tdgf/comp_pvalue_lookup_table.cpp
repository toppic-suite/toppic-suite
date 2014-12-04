#include <boost/algorithm/string.hpp>

#include "base/logger.hpp"
#include "tdgf/comp_pvalue_lookup_table.hpp"
#include "base/file_util.hpp"

namespace prot {

CompPValueLookupTable::CompPValueLookupTable(TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  initTable();

  test_num_ptr_ = CountTestNumPtr(new CountTestNum(mng_ptr));
  LOG_DEBUG("test number initialized");
}

void CompPValueLookupTable::initTable() {
  // add init table
  int ppo = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()
      ->getPpo();
  std::string line;
  std::vector<std::string> strs;

  input_.open(
      mng_ptr_->prsm_para_ptr_->getExeDir() + FILE_SEPARATOR
          + "toppic_resources" + FILE_SEPARATOR + "p_value_table"
          + FILE_SEPARATOR + "ppm" + std::to_string(ppo) + "_ptm0.txt",
      std::ios::in);

  while (std::getline(input_, line)) {
    boost::split(strs, line, boost::is_any_of(" \t"));
    ptm0_[getPeakIndex(std::stoi(strs[0]))][getFragIndex(std::stoi(strs[1]))] =
        std::stod(strs[2]);
  }

  input_.open(
      mng_ptr_->prsm_para_ptr_->getExeDir() + FILE_SEPARATOR
          + "toppic_resources" + FILE_SEPARATOR + "p_value_table"
          + FILE_SEPARATOR + "ppm" + std::to_string(ppo) + "_ptm1.txt",
      std::ios::in);

  while (std::getline(input_, line)) {
    boost::split(strs, line, boost::is_any_of(" \t"));
    ptm1_[getPeakIndex(std::stoi(strs[0]))][getFragIndex(std::stoi(strs[1]))] =
        std::stod(strs[2]);
  }

  input_.open(
      mng_ptr_->prsm_para_ptr_->getExeDir() + FILE_SEPARATOR
          + "toppic_resources" + FILE_SEPARATOR + "p_value_table"
          + FILE_SEPARATOR + "ppm" + std::to_string(ppo) + "_ptm2.txt",
      std::ios::in);

  while (std::getline(input_, line)) {
    boost::split(strs, line, boost::is_any_of(" \t"));
    ptm2_[getPeakIndex(std::stoi(strs[0]))][getFragIndex(std::stoi(strs[1]))] =
        std::stod(strs[2]);
  }

  LOG_DEBUG("table initialized");
}

double CompPValueLookupTable::compProb(int ppo, double prec_mass, int peak_num,
                                       int match_frag_num,
                                       int unexpected_shift_num) {
  // add implementation.
  return 1.0;
}

/* set alignment */
void CompPValueLookupTable::process(DeconvMsPtr deconv_ms_ptr,
                                    PrsmPtrVec &prsm_ptrs) {
  int ppo = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()
      ->getPpo();
  int peak_num = deconv_ms_ptr->size();
  double tolerance = deconv_ms_ptr->getHeaderPtr()->getErrorTolerance();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    double refine_prec_mass = prsm_ptrs[i]->getAdjustedPrecMass();
    int match_frag_num = prsm_ptrs[i]->getMatchFragNum();
    int unexpected_shift_num = prsm_ptrs[i]->getProteoformPtr()
        ->getUnexpectedChangeNum();

    double prot_prob = compProb(ppo, refine_prec_mass, peak_num, match_frag_num,
                                unexpected_shift_num);

    SemiAlignTypePtr type_ptr = prsm_ptrs[i]->getProteoformPtr()
        ->getSemiAlignType();

    double cand_num = test_num_ptr_->compCandNum(type_ptr, unexpected_shift_num,
                                                 refine_prec_mass, tolerance);

    ExtremeValuePtr prob_ptr(new ExtremeValue(prot_prob, cand_num, 1));

    prsm_ptrs[i]->setProbPtr(prob_ptr);

  }
}

int getPeakIndex(int p) {
  int k, t;
  if (p <= 10) {
    return 0;
  } else if (p <= 100) {
    k = p / 5;
    t = p % 5;
    if (t >= 3) {
      return k - 1;
    } else {
      return k - 2;
    }
  } else if (p <= 200) {
    p = p - 100;
    k = p / 10;
    t = p % 10;
    if (t >= 5) {
      return k + 19;
    } else {
      return k - 1 + 19;
    }
  } else if (p <= 400) {
    p = p - 200;
    k = p / 20;
    t = p % 20;
    if (t >= 10) {
      return k + 29;
    } else {
      return k - 1 + 29;
    }
  } else {
    if (p <= 425) {
      return 38;
    } else if (p <= 475)
      return 39;
    else
      return 40;
  }
}

int getFragIndex(int i) {
  if (i <= 5)
    return 0;

  int k, t;
  k = i / 5;
  t = i % 5;
  if (t >= 3)
    return k;
  else
    return k - 1;
}

}
