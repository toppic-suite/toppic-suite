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
  int ppo = mng_ptr_->prsm_para_ptr_->getErrorTolerance();
  LOG_DEBUG("ppo " << ppo);

  memset(ptm0_, 0, sizeof(ptm0_[0][0]) * 41 * 20);
  memset(ptm1_, 0, sizeof(ptm0_[0][0]) * 41 * 20);
  memset(ptm2_, 0, sizeof(ptm0_[0][0]) * 41 * 20);

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

  input_.close();

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

  input_.close();

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

  input_.close();
  LOG_DEBUG("table initialized");

}

double CompPValueLookupTable::compProb(int peak_num, int match_frag_num,
                                       int unexpected_shift_num) {
  // add implementation.
  std::vector<int> idx = getFourIndex(peak_num, match_frag_num);
  double res = 0;

  int x1 = idx[0], x2 = idx[1], y1 = idx[2], y2 = idx[3];

  double p11, p12, p21, p22;

  if (unexpected_shift_num == 0) {
    p11 = ptm0_[x1][y1];
    p12 = ptm0_[x1][y2];
    p21 = ptm0_[x2][y1];
    p22 = ptm0_[x2][y2];
  } else if (unexpected_shift_num == 1) {
    p11 = ptm1_[x1][y1];
    p12 = ptm1_[x1][y2];
    p21 = ptm1_[x2][y1];
    p22 = ptm1_[x2][y2];
  } else {
    p11 = ptm2_[x1][y1];
    p12 = ptm2_[x1][y2];
    p21 = ptm2_[x2][y1];
    p22 = ptm2_[x2][y2];
  }

  p11 = log(p11);
  p12 = log(p12);
  p21 = log(p21);
  p22 = log(p22);

  x1 = getPeakNumFromIndex(idx[0]);
  x2 = getPeakNumFromIndex(idx[1]);
  y1 = 5 * (idx[2] + 1);
  y2 = 5 * (idx[3] + 1);

  res = ((x2 - peak_num) * (y2 - match_frag_num) * p11
      + (peak_num - x1) * (y2 - match_frag_num) * p21
      + (x2 - peak_num) * (match_frag_num - y1) * p12
      + (peak_num - x1) * (match_frag_num - y1) * p22)
      / ((x2 - x1) * (y2 - y1));

  res = exp(res);

  LOG_DEBUG("prob " << res);

  return res;
}

/* set alignment */
void CompPValueLookupTable::process(const DeconvMsPtrVec &deconv_ms_ptr_vec, PrsmPtrVec &prsm_ptrs) {
  int ppo = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
  int peak_num = 0; 
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    peak_num += deconv_ms_ptr_vec[i]->size();
  }
  double tolerance = deconv_ms_ptr_vec[0]->getHeaderPtr()->getErrorTolerance();
  double refine_prec_mass = deconv_ms_ptr_vec[0]->getHeaderPtr()->getPrecMonoMassMinusWater();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    
    int match_frag_num = prsm_ptrs[i]->getMatchFragNum();
    int unexpected_shift_num = prsm_ptrs[i]->getProteoformPtr()
        ->getUnexpectedChangeNum();

    double prot_prob = 1.0;

    if (match_frag_num < 5) {
      prot_prob = 1.0;
    } else {
      prot_prob = compProb(peak_num, match_frag_num, unexpected_shift_num);
    }

    SemiAlignTypePtr type_ptr = prsm_ptrs[i]->getProteoformPtr()
        ->getSemiAlignType();

    double cand_num = test_num_ptr_->compCandNum(type_ptr, unexpected_shift_num,
                                                 refine_prec_mass, tolerance);

    ExtremeValuePtr prob_ptr(new ExtremeValue(prot_prob, cand_num, 1));

    prsm_ptrs[i]->setProbPtr(prob_ptr);

  }
}

bool CompPValueLookupTable::inTable(DeconvMsPtr deconv_ms_ptr,
                                    PrsmPtrVec &prsm_ptrs) {

  int peak_num = deconv_ms_ptr->size();

  if (peak_num > 500 || peak_num < 10)
    return false;

  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int match_frag_num = prsm_ptrs[i]->getMatchFragNum();

    if (match_frag_num > 100)
      return false;

    std::vector<int> idx = getFourIndex(peak_num, match_frag_num);

    int unexpected_shift_num = prsm_ptrs[i]->getProteoformPtr()
        ->getUnexpectedChangeNum();

    if (unexpected_shift_num == 0) {
      if (ptm0_[idx[0]][idx[2]] == 0 || ptm0_[idx[0]][idx[3]] == 0
          || ptm0_[idx[1]][idx[2]] == 0 || ptm0_[idx[1]][idx[3]] == 0)
        return false;
    } else if (unexpected_shift_num == 1) {
      if (ptm1_[idx[0]][idx[2]] == 0 || ptm1_[idx[0]][idx[3]] == 0
          || ptm1_[idx[1]][idx[2]] == 0 || ptm1_[idx[1]][idx[3]] == 0)
        return false;
    } else {
      if (ptm2_[idx[0]][idx[2]] == 0 || ptm2_[idx[0]][idx[3]] == 0
          || ptm2_[idx[1]][idx[2]] == 0 || ptm2_[idx[1]][idx[3]] == 0)
        return false;
    }
  }

  return true;

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

std::vector<int> getFourIndex(int peak_num, int frag_num) {
  std::vector<int> idx(4);

  int k;

  if (peak_num <= 10) {
    idx[0] = 0;
    idx[1] = 1;
  } else if (peak_num <= 100) {
    k = peak_num / 5;

    idx[0] = k - 2;
    idx[1] = k - 1;
  } else if (peak_num <= 200) {
    peak_num = peak_num - 100;
    k = peak_num / 10;
    idx[0] = k + 18;
    idx[1] = k + 19;
  } else if (peak_num <= 400) {
    peak_num = peak_num - 200;
    k = peak_num / 20;
    idx[0] = k + 28;
    idx[1] = k + 29;
  } else {
    if (peak_num <= 450) {
      idx[0] = 38;
      idx[1] = 39;
    } else {
      idx[0] = 39;
      idx[1] = 40;
    }
  }

  if (frag_num <= 5) {
    idx[2] = 0;
    idx[3] = 1;
  } else {

    k = frag_num / 5;

    idx[2] = k - 1;
    idx[3] = k;
  }
  return idx;
}

int getPeakNumFromIndex(int idx) {
  if (idx <= 18) {
    return 5 * (idx + 2);
  } else if (idx <= 28) {
    return 100 + 10 * (idx - 18);
  } else if (idx <= 38) {
    return 200 + 20 * (idx - 28);
  } else {
    return 400 + 50 * (idx - 38);
  }
}

}
