#ifndef PROT_LOCAL_UTIL_HPP_
#define PROT_LOCAL_UTIL_HPP_

#include "base/ptm.hpp"
#include "base/algorithm.hpp"
#include "prsm/prsm_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "spec/theo_peak.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "local_mng.hpp"

namespace prot {

void getNtermTruncRange(PrsmPtr & prsm, int & min, int & max, double max_mass);

void getCtermTruncRange(PrsmPtr & prsm, int & min, int & max, double max_mass);

void scr_filter(const std::vector<double> & scr, int & bgn, int & end,
                double & conf, double thread);

void getSupPeakNum(const PrsmPtr & prsm, const ChangePtr & change,
                   double min_mass, int & left, int &right);

bool modifiable(const ProteoformPtr& proteoform_ptr, int i,
                const PtmPtr& ptm_ptr, bool cysteine_protected);

template<typename T>
std::vector<double> normalize(const std::vector<T>& scr,
                              typename std::enable_if<std::is_arithmetic<T>::value>::type* = 0) {
    std::vector<double> res(scr.size());
    double max = *std::max_element(scr.begin(), scr.end());
    for (size_t i = 0; i < scr.size(); i++) {
        res[i] = scr[i] / max;
    }

    double sum = std::accumulate(res.begin(), res.end(), 0.0);
    for (size_t i = 0; i < res.size(); i++) {
        res[i] /= sum;
    }
    return res;
}

void termAdjust(PrsmPtr& prsm, bool& ptm_known, double& t_one_ptm_scr, const LocalMngPtr& mng_ptr);

void termAdjust(PrsmPtr& prsm, bool& ptm1_known, bool &ptm2_known, double& t_two_ptm_scr, const LocalMngPtr& mng_ptr);

void massAdjust(PrsmPtr& prsm, const LocalMngPtr& mng_ptr);

double getScr(PrsmPtr& prsm, const LocalMngPtr& mng_ptr);

double getScr(PrsmPtr& prsm, bool known, const LocalMngPtr& mng_ptr);

int getSplit(const PrsmPtr& prsm, const PtmPtr& ptm1, const PtmPtr& ptm2, const LocalMngPtr& mng_ptr);

int getNumMatch(double theo_mod, const std::vector<double>& spec_double,
                const std::vector<double>& spec_torlerance);

} //namespace prot
#endif /* PROT_LOCAL_UTIL_HPP */
