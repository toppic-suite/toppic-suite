/*
 * util.cpp
 *
 *  Created on: May 14, 2014
 *      Author: kouqiang
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <map>
#include <cstring>
#include <algorithm>

#include "local/util.hpp"
#include "base/ptm.hpp"

namespace prot {

int Util::ppm = 15;
double Util::preMass = 12000.0;
double Util::cysteineProtection = 0.0;

Util::Util(int a, double preMass, std::string protection) {
  this->ppm = a;
  this->preMass = preMass;
  if (protection == "C0") {
    cysteineProtection = 0.0;
  } else if (protection == "C57") {
    cysteineProtection = 57.021464;
  } else {
    cysteineProtection = 58.005479;
  }
}

Util::Util(double preMass) {
  this->preMass = preMass;
}

bool Util::InErrorRange(double a, double b) {
  double thread = preMass * ppm / 1000000;
//  double thread = b * ppm / 1000000;
  if (thread <= 0.01)
    thread = 0.01;

  if (std::abs(a - b) <= thread)
    return true;
  else if (std::abs(std::abs(a - b) - 1) <= thread)
    return true;
  else
    return false;
}

/*
 * get the parameters from the PrSM without PTM
 */
std::vector<double> Util::getProb(std::vector<double> theo, PrmPeakPtrVec prm) {
  std::vector<double> res(4, 0.0);
  double t1, t3 = 0;
  bool match = false;
  for (size_t i = 0; i < prm.size(); i++) {
    match = false;
    for (size_t j = 0; j < theo.size(); j++) {
      t1 = theo[j];
      if (std::abs(theo[j] - prm[i]->getMonoMass())
          < prm[i]->getStrictTolerance()) {
        res[2] += (theo[j] * 2 * ppm / 1000000);
        match = true;
        break;
      }
    }
    if (!match)
      res[0] += t1 * 2 * ppm / 1000000;

  }

  for (size_t i = 0; i < prm.size(); i++) {
    t3 += prm[i]->getMonoMass() * 2 * ppm / 1000000;
  }
  res[1] = t3 - res[2];
  res[3] = preMass;

  return res;
}

std::vector<double> Util::getPara(DeconvMsPtrVec msalign, PrsmPtrVec prsms,
                                  SpParaPtr sp_para) {
  std::vector<double> para, res(4, 0.0);
  std::vector<double> theo;
  std::string pep;
  DeconvMsPtr DeconvRes;
  for (size_t i = 0; i < prsms.size(); i++) {
    if (prsms[i]->getProteoformPtr()->getChangePtrVec().size() == 0) {
      DeconvRes = getDeconvMsPtrbyID(msalign, prsms[i]->getSpectrumId());
      SpectrumSetPtr spec_set_ptr = getSpectrumSet(DeconvRes, 0.0, sp_para);
      PrmPeakPtrVec prm = spec_set_ptr->getSpSix()->getPeakPtrVec();
      pep = prsms[i]->getProteoformPtr()->getResSeqPtr()->toString();
      theo = getTheo(pep.substr(0, pep.length() - 1));
      para = getProb(theo, prm);
      for (int j = 0; j < 4; j++) {
        res[j] += para[j];
      }
    }
  }
  return res;
}

double Util::getPTMProb(std::vector<double> theo, PrmPeakPtrVec prm, double p1,
                        double p2) {
  if (theo.size() == 0)
    return 0.0;

  double res;
  int t1 = theo.size(), t3 = 0;

  for (size_t i = 0; i < prm.size(); i++) {

    for (size_t j = 0; j < theo.size(); j++) {

      if (std::abs(prm[i]->getMonoMass() - theo[j])
          < prm[i]->getStrictTolerance()) {
        t3++;
      }
    }
  }
  t1 = theo.size();
  res = pow(p1, t1) * pow(p2, t3);
  return res;

}
/*
 * to modify a theoretical spectrum with a PTM
 * k starts from 0
 */
std::vector<double> Util::modSeq(std::vector<double> theo, int k, double mod) {
  int n = round(abs(mod) / 110.0 + 1);
  if (mod + 57.0 > 0) {
    for (size_t i = k; i < theo.size(); i++) {
      theo[i] = theo[i] + mod;
    }
  } else {
    if (n + k >= (int) theo.size()) {
      return std::vector<double>();
    } else {
      for (size_t i = k; i < theo.size(); i++) {
        theo[i] = theo[i] + mod;
      }
    }
  }

  for (size_t i = 0; i < theo.size(); i++) {
    if (theo[i] <= 50.0)
      theo[i] = 0.0;
  }

  return theo;
}

/*
 * get a theoretical spectrum from a peptide only the prefix
 */
std::vector<double> Util::getTheo(std::string peptide) {
  std::vector<double> theo;
  for (size_t i = 1; i <= peptide.length(); i++) {
    std::string begin = peptide.substr(0, i);
    theo.push_back(AcidFactory::getPeptideMass(begin, cysteineProtection));
  }
  std::sort(theo.begin(), theo.end());
  return theo;
}
// k starts from 1
std::vector<double> Util::getTheo(std::string peptide, int k) {
  std::vector<double> theo;
  for (int i = 1; i <= k; i++) {
    std::string begin = peptide.substr(0, i);
    theo.push_back(AcidFactory::getPeptideMass(begin, cysteineProtection));
  }
  std::sort(theo.begin(), theo.end());
  return theo;
}

std::vector<double> Util::getTheo(std::string peptide, int j, int k) {
  std::vector<double> theo;
  for (int i = j; i <= k; i++) {
    std::string begin = peptide.substr(0, i);
    theo.push_back(AcidFactory::getPeptideMass(begin, cysteineProtection));
  }
  std::sort(theo.begin(), theo.end());
  return theo;
}

std::vector<double> Util::getTheoByPos(std::string peptide, int pos,
                                       std::vector<double> theo) {
  // pos start from 1, not 0

  std::string begin = peptide.substr(0, pos);
  theo.push_back(AcidFactory::getPeptideMass(begin));

  std::string end = peptide.substr(pos, peptide.length() - pos);
  theo.push_back(AcidFactory::getPeptideMass(end) + 18);

  return theo;
}

std::vector<double> Util::getTheoByPos(std::string peptide,
                                       std::vector<double> theo) {
  // pos start from 1, not 0

  std::string begin = peptide.substr(0, 1);
  theo.push_back(AcidFactory::getPeptideMass(begin));
  theo.push_back(AcidFactory::getPeptideMass(peptide));
  std::string end = peptide.substr(1, peptide.length() - 1);
  theo.push_back(AcidFactory::getPeptideMass(end) + 18);

  return theo;
}

/*
 * all possible combination of PTM by two mass value
 * return a vector
 */

PtmPtrVec Util::getPTM(double a, double b) {
  PtmPtrVec res;
  double mass = a + b;
  PtmPtrCom ptms = PtmFactory::getPtmPtrCombine();
  for (size_t i = 0; i < ptms.size(); i++) {
    if (InErrorRange(mass, ptms[i].first)
        && fabs(ptms[i].second.first->getMonoMass()) > 1
        && fabs(ptms[i].second.second->getMonoMass()) > 1) {
      res.push_back(ptms[i].second.first);
      res.push_back(ptms[i].second.second);
    }
  }
  return res;
}

PtmPtrVec Util::getPTM(double mass) {
  PtmPtrVec res;

  PtmPtrVec ptms = PtmFactory::getBasePtmPtrVec();
  for (size_t i = 0; i < ptms.size(); i++) {
    if (InErrorRange(mass, ptms[i]->getMonoMass())
        && fabs(ptms[i]->getMonoMass()) > 1) {
      res.push_back(ptms[i]);
    }
  }
  return res;
}

PtmPtrVec Util::getPTM(double a, double b, double c) {
  PtmPtrVec res;
  double mass = a + b + c;
  PtmPtrCom ptms = PtmFactory::getPtmPtrCombine();
  for (size_t i = 0; i < ptms.size(); i++) {
    if (InErrorRange(mass, ptms[i].first)
        && fabs(ptms[i].second.first->getMonoMass()) > 1
        && fabs(ptms[i].second.second->getMonoMass()) > 1) {
      res.push_back(ptms[i].second.first);
      res.push_back(ptms[i].second.second);
    }
  }
  return res;
}

int Util::getMatch(double peak, std::vector<double> exp) {
  int res = 0;
  double theo;
  for (size_t i = 0; i < exp.size(); i++) {
    theo = exp[i];
    if (InErrorRange(theo, peak)) {
      res++;
    }
  }
  return res;
}

int Util::getMatch(std::vector<double> theo, std::vector<double> exp) {
  if (theo.size() == 0)
    return 0;
  else {
    double expMass, theoMass;
    int t3 = 0;
    std::vector<double>::iterator i, j;
    for (i = exp.begin(); i != exp.end(); i++) {
      expMass = *i;
      for (j = theo.begin(); j != theo.end(); j++) {
        theoMass = *j;
        if (InErrorRange(expMass, theoMass)) {
          t3++;
        }
      }
    }
    return t3;
  }

}

//double Util::addPartTheo(std::string peptide, int k) {
//	std::string begin = peptide.substr(0, k);
//	return AcidFactory::getPeptideMass(begin);
//}

std::vector<int> Util::getScoreVec(std::vector<double> scr, double a) {
  std::vector<int> res;

  return res;
}

bool Util::changeFilter(ProteoformPtr proteoform, int i) {
  bool res = true;
  ChangePtr change = proteoform->getChangePtrVec()[i];
  if (change->getMassShift() == cysteineProtection) {
    res = false;
  } else if (abs(change->getMassShift()) <= 1.0) {
    res = false;
  } else if (change->getChangeType() == 1 || change->getChangeType() == 2) {
    res = false;
  }
  return res;
}

int Util::changeFilter(ProteoformPtr proteoform) {
  int count = 0;
  ChangePtrVec change = proteoform->getChangePtrVec();
  for (size_t i = 0; i < change.size(); i++) {
    if (change[i]->getChangeType() != 1 && change[i]->getChangeType() != 2
        && abs(change[i]->getMassShift()) > 1
        && change[i]->getMassShift() != cysteineProtection)
      count++;
  }
  return count;
}

std::vector<int> Util::changeFilterID(ProteoformPtr proteoform) {
  std::vector<int> res;
  ChangePtrVec change = proteoform->getChangePtrVec();
  for (size_t i = 0; i < change.size(); i++) {
    if (change[i]->getChangeType() != 1 && change[i]->getChangeType() != 2
        && abs(change[i]->getMassShift()) > 1
        && change[i]->getMassShift() != cysteineProtection)
      res.push_back(i);
  }
  return res;
}

std::string Util::posFilter(std::string pos) {
  std::string res;
  if (cysteineProtection > 0) {
    for (size_t i = 0; i < pos.size(); ++i) {
      if (pos[i] != 'C')
        res += pos[i];
    }
    return res;
  } else
    return pos;

}

}

