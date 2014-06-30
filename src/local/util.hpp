/*
 * util.hpp
 *
 *  Created on: May 14, 2014
 *      Author: kouqiang
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include "base/base_data.hpp"
#include "base/ptm.hpp"
#include "base/acid.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_combine.hpp"
#include "prsm/prsm_selector.hpp"
#include "prsm/output_selector.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/table_writer.hpp"
#include "prsm/prsm_fdr.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/prm_peak.hpp"
#include "spec/spectrum_set.hpp"

#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace prot {

class Util {
 public:
  Util(int ppm, double preMass, std::string protection);

  Util(double preMass);

  static bool InErrorRange(double a, double b);

  static std::vector<double> getProb(std::vector<double> theo,
                                     PrmPeakPtrVec prm);

  static double getPTMProb(std::vector<double> theo, PrmPeakPtrVec prm,
                           double p1, double p2);

  static double getPTMProb(std::vector<double> theo, std::vector<double> spec,
                           double p1, double p2, int t);

  static std::vector<double> modSeq(std::vector<double> theo, int k,
                                    double mod);

  static std::vector<double> getTheo(std::string peptide);

  static std::vector<double> getTheo(std::string peptide, int k);

  static std::vector<double> getTheo(std::string peptide, int j, int k);

  static std::vector<double> getTheoByPos(std::string peptide, int pos,
                                          std::vector<double> theo);

  static std::vector<double> getTheoByPos(std::string peptide,
                                          std::vector<double> theo);

  static PtmPtrVec getPTM(double mass);

  static PtmPtrVec getPTM(double a, double b);

  static PtmPtrVec getPTM(double a, double b, double c);

  static std::vector<double> getPartTheo(std::string peptide, int l);

  static int getMatch(double peak, std::vector<double> spec);

  static int getMatch(std::vector<double> theo, std::vector<double> spec);

  static std::vector<int> getScoreVec(std::vector<double> scr, double a);

  static std::vector<double> getPara(DeconvMsPtrVec msalign, PrsmPtrVec prsms,
                                     SpParaPtr sp_para);

  static bool changeFilter(ProteoformPtr proteoform, int i);

  static int changeFilter(ProteoformPtr proteoform);

  static std::string posFilter(std::string pos);

  static std::vector<int> changeFilterID(ProteoformPtr proteoform);

 private:

  static int ppm;
  static double preMass;
  static double cysteineProtection;
};

inline std::vector<double> comScore(std::vector<double> s1,
                                    std::vector<double> s2, double w1,
                                    double w2) {
  if (s1.size() != s2.size()) {
    std::cerr << "Input vectors have different length!" << std::endl;
    exit(1);
  } else {
    std::vector<double> res;
    for (size_t i = 0; i < s1.size(); i++) {
      res.push_back(s1[i] * w1 + s2[i] * w2);
    }
    return res;
  }
}

inline std::vector<double> getMaxVec(std::vector<std::vector<double> > scrVec) {
  std::vector<double> res, sum;
  double temp;
  int index;
  for (size_t i = 0; i < scrVec.size(); i++) {
    temp = std::accumulate(scrVec[i].begin(), scrVec[i].end(), 0.0);
    sum.push_back(temp);
  }
  index = std::distance(sum.begin(), std::max_element(sum.begin(), sum.end()));
  return scrVec[index];
}

inline int getMaxIndex(std::vector<std::vector<double> > scrVec) {
  std::vector<double> sum;
  double temp;
  int index;
  for (size_t i = 0; i < scrVec.size(); i++) {
    temp = std::accumulate(scrVec[i].begin(), scrVec[i].end(), 0.0);
    sum.push_back(temp);
  }
  index = std::distance(sum.begin(), std::max_element(sum.begin(), sum.end()));
  return index;
}

inline int getBreakPoint(std::vector<std::vector<double> > mat) {
  int res;
  std::vector<double> sumVec;
  double sum = 0.0;
  for (size_t i = 1; i < mat.size(); i++) {
    for (size_t j = 0; j < mat.size(); j++) {
      for (size_t k = 0; k < mat.size(); k++) {
        if (j <= i && k >= i) {
          sum += mat[j][k];
        }
      }
    }
    sumVec.push_back(sum);
    sum = 0;
  }
  res = std::distance(sumVec.begin(),
                      std::max_element(sumVec.begin(), sumVec.end()));
  return res;
}

inline double getPosProb(std::vector<std::vector<double> > mat) {
  double res;
  std::vector<double> sumVec;
  double sum = 0.0;
  for (size_t i = 1; i < mat.size(); i++) {
    for (size_t j = 0; j < mat.size(); j++) {
      for (size_t k = 0; k < mat.size(); k++) {
        if (j <= i && k >= i) {
          sum += mat[j][k];
        }
      }
    }
    sumVec.push_back(sum);
    sum = 0;
  }
  res = *std::max_element(sumVec.begin(), sumVec.end())
      / std::accumulate(sumVec.begin(), sumVec.end(), 0.0);
  return res;
}

inline double getMassProb(std::vector<std::vector<double> > scrVec) {
  std::vector<double> sum;
  double temp;
  for (size_t i = 0; i < scrVec.size(); i++) {
    temp = std::accumulate(scrVec[i].begin(), scrVec[i].end(), 0.0);
    sum.push_back(temp);
  }
  temp = *std::max_element(sum.begin(), sum.end());
  return temp / (std::accumulate(sum.begin(), sum.end(), 0.0));
}

inline std::vector<double> scaleVec(std::vector<double> scr) {
  std::vector<double> res1, res2;
  double max = *std::max_element(scr.begin(), scr.end());
  for (size_t i = 0; i < scr.size(); i++) {
    res1.push_back(scr[i] / max);
  }

  double sum = std::accumulate(res1.begin(), res1.end(), 0.0);

  for (size_t i = 0; i < res1.size(); i++) {
    res2.push_back(res1[i] / sum);
  }
  return res2;
}

inline DeconvMsPtr getDeconvMsPtrbyID(DeconvMsPtrVec msalign, int id) {
  DeconvMsPtr res(nullptr);
  for (size_t i = 0; i < msalign.size(); i++) {
    if (msalign[i]->getHeaderPtr()->getId() == id)
      res = msalign[i];
  }
  return res;
}

}

#endif /* UTIL_HPP_ */
