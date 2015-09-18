#include <boost/range/algorithm/remove_if.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "spec/spectrum_set.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "local_util.hpp"

namespace prot {

void scr_filter(const std::vector<double> & scr, int & bgn, int & end,
                double & conf, double thread) {

    bgn = std::distance(scr.begin(), std::max_element(scr.begin(), scr.end()));
    end = scr.size();

    for (; end >= 0; end--) {
        if (scr[bgn] == scr[end]) {
            break;
        }
    }

    conf = 0;

    for (int i = bgn; i <= end; i++) {
        conf += scr[i];
    }

    while (conf <= thread) {
        if (bgn == 0) {
            end++;
            conf += scr[end];
        } else if (end == (int) scr.size() - 1) {
            bgn--;
            conf += scr[bgn];
        } else {
            if (scr[bgn] == scr[end]) {
                end++;
                conf += scr[end];
                bgn--;
                conf += scr[bgn];
            } else if (scr[bgn] > scr[end]) {
                bgn--;
                conf += scr[bgn];
            } else {
                end++;
                conf += scr[end];
            }
        }
    }
}

void getNtermTruncRange(PrsmPtr & prsm, int & min, int & max, double max_mass) {
    double ori_mass = prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0]
                      ->getMassShift();
    double mass = ori_mass;

    int ori_start = prsm->getProteoformPtr()->getStartPos();

    max = min = 0;

    while (ori_start + min > 0) {
        min--;
        std::string t_seq =
            prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                ori_start + min, -min);
        mass = ori_mass - AcidFactory::getPeptideMass(t_seq);

        if (std::abs(mass) >= max_mass || ori_start + min <= 0) {
            min++;
            break;
        }
    }

    mass = ori_mass;
    while (std::abs(ori_mass) < max_mass) {
        max++;
        std::string t_seq =
            prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                ori_start, max);
        mass = ori_mass + AcidFactory::getPeptideMass(t_seq);

        if (std::abs(mass) >= max_mass) {
            max--;
            break;
        }
    }
}

void getCtermTruncRange(PrsmPtr & prsm, int & min, int & max, double max_mass) {

    double ori_mass = prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[prsm
                      ->getProteoformPtr()->getUnexpectedChangeNum() - 1]->getMassShift();
    double mass = ori_mass;

    int ori_end = prsm->getProteoformPtr()->getEndPos();

    max = min = 0;

    while (ori_end + std::abs(max)
            < prsm->getProteoformPtr()->getDbResSeqPtr()->getLen() - 1) {
        max++;
        std::string t_seq =
            prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                ori_end + 1, max);
        mass = ori_mass - AcidFactory::getPeptideMass(t_seq);

        if (std::abs(mass) >= max_mass
                || ori_end + std::abs(max)
                >= prsm->getProteoformPtr()->getDbResSeqPtr()->getLen() - 1) {
            max--;
            break;
        }

    }

    mass = ori_mass;

    while (std::abs(ori_mass) < max_mass) {
        min--;
        std::string t_seq =
            prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                ori_end + 1 + min, -min);
        mass = ori_mass + AcidFactory::getPeptideMass(t_seq);
        if (std::abs(mass) >= max_mass) {
            min++;
            break;
        }
    }
}

void getSupPeakNum(const PrsmPtr & prsm, const ChangePtr & change, double min_mass, int & left, int &right) {
    left = right = 0;
    PeakIonPairPtrVec pair_ptrs = getPeakIonPairs(prsm->getProteoformPtr(),
                                  prsm->getRefineMsPtrVec(), min_mass);

    for (size_t i = 0; i < pair_ptrs.size(); i++) {
        if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getIonTypePtr()->getName() == "B"
                || pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getIonTypePtr()->getName() == "C") {
            if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getDisplayPos()
                    <= change->getLeftBpPos()) {
                left++;
            }
            if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getDisplayPos()
                    >= change->getRightBpPos()) {
                right++;
            }
        } else {
            if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getDisplayPos()
                    >= prsm->getProteoformPtr()->getLen() - change->getLeftBpPos()) {
                left++;
            }
            if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getDisplayPos()
                    <= prsm->getProteoformPtr()->getLen() - change->getRightBpPos()) {
                right++;
            }
        }
    }
}

bool modifiable(const ProteoformPtr& proteoform_ptr, int i,
                const PtmPtr& ptm_ptr, bool cysteine_protected) {

    for (size_t j = 0; j < proteoform_ptr->getChangePtrVec().size(); j++) {
        if (proteoform_ptr->getChangePtrVec()[j]->getChangeType() != UNEXPECTED_CHANGE) {
            if (i < proteoform_ptr->getChangePtrVec()[j]->getRightBpPos()
                    && i >= proteoform_ptr->getChangePtrVec()[j]->getLeftBpPos()) {
                return false;
            }
        }
    }

    if (ptm_ptr == nullptr) {
        return true;
    }

    int start = proteoform_ptr->getStartPos();
    int end = proteoform_ptr->getEndPos();

    std::string posN = ptm_ptr->getPosN();
    std::string posC = ptm_ptr->getPosC();
    std::string pos = ptm_ptr->getPos();

    if (cysteine_protected) {
        posN.erase(boost::remove_if(posN, boost::is_any_of("C")), posN.end());
        posC.erase(boost::remove_if(posC, boost::is_any_of("C")), posC.end());
        pos.erase(boost::remove_if(pos, boost::is_any_of("C")), pos.end());
    }

    std::string pep = proteoform_ptr->getResSeqPtr()->toAcidString();

    std::size_t found;

    if (i == 0 && start == 0) {
        std::string temp = posN + pos;
        found = temp.find(std::string(1, pep[i]));
    } else if (i + start == end
               && end == (proteoform_ptr->getDbResSeqPtr()->getLen() - 1)) {
        std::string temp = posC + pos;
        found = temp.find(std::string(1, pep[i]));
    } else {
        found = pos.find(std::string(1, pep[i]));
    }

    return found != std::string::npos;

}

void termAdjust(PrsmPtr& prsm, bool& ptm_known, double& t_one_ptm_scr, const LocalMngPtr& mng_ptr) {
    double ppm = mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
    int left_sup, right_sup;
    ChangePtr change_ptr =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
    int ori_start = prsm->getProteoformPtr()->getStartPos();
    int ori_end = prsm->getProteoformPtr()->getEndPos();

    double ori_mass = change_ptr->getMassShift();
    double mass = ori_mass;
    getSupPeakNum(prsm, change_ptr, mng_ptr->min_mass_, left_sup, right_sup);

    LOG_DEBUG("left_sup " << left_sup << " right_sup " << right_sup);

    int n_trunc_min, n_trunc_max, c_trunc_min, c_trunc_max;
    getNtermTruncRange(prsm, n_trunc_min, n_trunc_max, mng_ptr->max_ptm_mass_);
    getCtermTruncRange(prsm, c_trunc_min, c_trunc_max, mng_ptr->max_ptm_mass_);

    LOG_DEBUG("n_trunc_min " << n_trunc_min << " n_trunc_max " << n_trunc_max);
    LOG_DEBUG("c_trunc_min " << n_trunc_min << " c_trunc_max " << n_trunc_max);

    std::vector<bool> ptm_known_vec;
    std::vector<double> t_scr;
    std::vector<int> c_vec;
    std::vector<int> n_vec;
    std::string n_seq, c_seq;

    int n_t = 0, n_f = 0;

    if (left_sup <= LEFT_SUP_LIMIT && right_sup <= RIGHT_SUP_LIMIT) {
        for (int i = n_trunc_min; i <= n_trunc_max; i++) {
            for (int j = c_trunc_min; j <= c_trunc_max; j++) {
                if (i < 0) {
                    n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                            .substr(ori_start + i, -i);
                    mass = ori_mass - AcidFactory::getPeptideMass(n_seq);
                } else {
                    n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                            .substr(ori_start, i);
                    mass = ori_mass + AcidFactory::getPeptideMass(n_seq);
                }

                if (j >= 0) {
                    c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                            .substr(ori_end + 1, j);
                    mass = mass - AcidFactory::getPeptideMass(c_seq);
                } else {
                    c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                            .substr(ori_end + 1 + j, -j);
                    mass = mass + AcidFactory::getPeptideMass(c_seq);
                }

                if (std::abs(mass) > mng_ptr->max_ptm_mass_) {
                    continue;
                }

                n_vec.push_back(i);
                c_vec.push_back(j);

                change_ptr->setMassShift(mass);
                prsm->getProteoformPtr()->setStartPos(ori_start + i);
                prsm->getProteoformPtr()->setEndPos(ori_end + j);

                if (std::abs(mass) <= 1 + prsm->getOriPrecMass() * ppm) {
                    return;
                }

                bool is_known = PtmFactory::isKnown(mass,
                                                    prsm->getOriPrecMass() * ppm);
                ptm_known_vec.push_back(is_known);

                t_scr.push_back(getScr(prsm, mng_ptr));
            }
        }

        if (ptm_known_vec.size() == 0) {
            ptm_known = false;
            return;
        }

        for (size_t i = 0; i < ptm_known_vec.size(); i++) {
            if (ptm_known_vec[i])
                n_t++;
            else
                n_f++;
        }

        for (size_t i = 0; i < t_scr.size(); i++) {
            if (ptm_known_vec[i]) {
                t_scr[i] = t_scr[i] * mng_ptr->theta_ / n_t;
            } else {
                t_scr[i] = t_scr[i] * (1 - mng_ptr->theta_) / n_f;
            }
        }

        int idx = std::distance(t_scr.begin(),
                                std::max_element(t_scr.begin(), t_scr.end()));
        ptm_known = ptm_known_vec[idx];
        t_one_ptm_scr = t_scr[idx];
        prsm->getProteoformPtr()->setStartPos(ori_start + n_vec[idx]);
        prsm->getProteoformPtr()->setEndPos(ori_end + c_vec[idx]);

        if (n_vec[idx] < 0) {
            n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_start + n_vec[idx], -n_vec[idx]);
            mass = ori_mass - AcidFactory::getPeptideMass(n_seq);
        } else {
            n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_start, n_vec[idx]);
            mass = ori_mass + AcidFactory::getPeptideMass(n_seq);
        }
        if (c_vec[idx] >= 0) {
            c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_end + 1, c_vec[idx]);
            mass = mass - AcidFactory::getPeptideMass(c_seq);
        } else {
            c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_end + 1 + c_vec[idx], -c_vec[idx]);
            mass = mass + AcidFactory::getPeptideMass(c_seq);
        }
        change_ptr->setMassShift(mass);

    } else if (left_sup <= LEFT_SUP_LIMIT) {

        for (int i = n_trunc_min; i <= n_trunc_max; i++) {
            n_vec.push_back(i);
            if (i < 0) {
                n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                        .substr(ori_start + i, -i);
                mass = ori_mass - AcidFactory::getPeptideMass(n_seq);
            } else {
                n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                        .substr(ori_start, i);
                mass = ori_mass + AcidFactory::getPeptideMass(n_seq);
            }

            change_ptr->setMassShift(mass);
            prsm->getProteoformPtr()->setStartPos(ori_start + i);

            if (std::abs(mass) <= 1 + prsm->getOriPrecMass() * ppm) {
                return;
            }

            bool is_known = PtmFactory::isKnown(mass, prsm->getOriPrecMass() * ppm);
            ptm_known_vec.push_back(is_known);
            t_scr.push_back(getScr(prsm, mng_ptr));
        }

        for (size_t i = 0; i < ptm_known_vec.size(); i++) {
            if (ptm_known_vec[i])
                n_t++;
            else
                n_f++;
        }

        for (size_t i = 0; i < t_scr.size(); i++) {
            if (ptm_known_vec[i]) {
                t_scr[i] = t_scr[i] * mng_ptr->theta_ / n_t;
            } else {
                t_scr[i] = t_scr[i] * (1 - mng_ptr->theta_) / n_f;
            }
        }

        int idx = std::distance(t_scr.begin(),
                                std::max_element(t_scr.begin(), t_scr.end()));
        ptm_known = ptm_known_vec[idx];
        t_one_ptm_scr = t_scr[idx];
        prsm->getProteoformPtr()->setStartPos(ori_start + n_vec[idx]);

        if (n_vec[idx] < 0) {
            n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_start + n_vec[idx], -n_vec[idx]);
            mass = ori_mass - AcidFactory::getPeptideMass(n_seq);
        } else {
            n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_start, n_vec[idx]);
            mass = ori_mass + AcidFactory::getPeptideMass(n_seq);
        }

        change_ptr->setMassShift(mass);

    } else if (right_sup <= RIGHT_SUP_LIMIT) {
        for (int i = c_trunc_min; i <= c_trunc_max; i++) {

            c_vec.push_back(i);

            if (i >= 0) {
                c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                        .substr(ori_end + 1, i);
                mass = ori_mass - AcidFactory::getPeptideMass(c_seq);
            } else {
                c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                        .substr(ori_end + 1 + i, -i);
                mass = ori_mass + AcidFactory::getPeptideMass(c_seq);
            }

            change_ptr->setMassShift(mass);
            prsm->getProteoformPtr()->setEndPos(ori_end + i);

            if (std::abs(mass) <= 1 + prsm->getOriPrecMass() * ppm) {
                return;
            }

            bool is_known = PtmFactory::isKnown(mass, prsm->getOriPrecMass() * ppm);
            ptm_known_vec.push_back(is_known);
            t_scr.push_back(getScr(prsm, mng_ptr));
        }

        for (size_t i = 0; i < ptm_known_vec.size(); i++) {
            if (ptm_known_vec[i])
                n_t++;
            else
                n_f++;
        }

        for (size_t i = 0; i < t_scr.size(); i++) {
            if (ptm_known_vec[i]) {
                t_scr[i] = t_scr[i] * mng_ptr->theta_ / n_t;
            } else {
                t_scr[i] = t_scr[i] * (1 - mng_ptr->theta_) / n_f;
            }
        }

        int idx = std::distance(t_scr.begin(),
                                std::max_element(t_scr.begin(), t_scr.end()));
        ptm_known = ptm_known_vec[idx];
        t_one_ptm_scr = t_scr[idx];
        prsm->getProteoformPtr()->setEndPos(ori_end + c_vec[idx]);

        if (c_vec[idx] >= 0) {
            c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_end + 1, c_vec[idx]);
            mass = ori_mass - AcidFactory::getPeptideMass(c_seq);
        } else {
            c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_end + 1 + c_vec[idx], -c_vec[idx]);
            mass = ori_mass + AcidFactory::getPeptideMass(c_seq);
        }
        change_ptr->setMassShift(mass);
    }
}

double getScr(PrsmPtr &prsm, const LocalMngPtr& mng_ptr) {

    if (prsm->getProteoformPtr()->getUnexpectedChangePtrVec().size() == 1) {

        double res = 0.0;
        ChangePtr change_ptr =
            prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
        int n = prsm->getProteoformPtr()->getLen();

        for (int j = 0; j < n; j++) {
            if (modifiable(prsm->getProteoformPtr(), j, nullptr, mng_ptr->cysteine_protected_)) {
                change_ptr->setLeftBpPos(j);
                change_ptr->setRightBpPos(j + 1);
                int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                               prsm->getRefineMsPtrVec(), mng_ptr->min_mass_);
                res += std::pow(mng_ptr->p1_, n - match) * std::pow(mng_ptr->p2_, match);
            }
        }

        return res;

    } else {

        double ppm = mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();

        ChangePtr change_ptr1 =
            prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
        ChangePtr change_ptr2 =
            prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[1];

        double mass1 = change_ptr1->getMassShift();
        double mass2 = change_ptr2->getMassShift();

        PtmPairVec ptm_pair_vec;

        bool ptm1_known = PtmFactory::isKnown(mass1, prsm->getOriPrecMass() * ppm);
        bool ptm2_known = PtmFactory::isKnown(mass2, prsm->getOriPrecMass() * ppm);

        if (ptm1_known && ptm2_known) {
            PtmPtrVec ptmVec1 = PtmFactory::getBasePtmPtrByMass(
                                    mass1, prsm->getOriPrecMass() * ppm);
            PtmPtrVec ptmVec2 = PtmFactory::getBasePtmPtrByMass(
                                    mass2, prsm->getOriPrecMass() * ppm);

            for (size_t i = 0; i < ptmVec1.size(); i++) {
                for (size_t j = 0; j < ptmVec2.size(); j++) {
                    ptm_pair_vec.push_back(std::make_pair(ptmVec1[i], ptmVec2[j]));
                }
            }

        } else if (ptm1_known) {
            PtmPtrVec ptmVec1 = PtmFactory::getBasePtmPtrByMass(
                                    mass1, prsm->getOriPrecMass() * ppm);
            for (size_t i = 0; i < ptmVec1.size(); i++) {
                ptm_pair_vec.push_back(std::make_pair(ptmVec1[i], nullptr));
            }
        } else if (ptm2_known) {
            PtmPtrVec ptmVec2 = PtmFactory::getBasePtmPtrByMass(
                                    mass2, prsm->getOriPrecMass() * ppm);
            for (size_t i = 0; i < ptmVec2.size(); i++) {
                ptm_pair_vec.push_back(std::make_pair(nullptr, ptmVec2[i]));
            }
        } else {
            ptm_pair_vec.push_back(std::make_pair(nullptr, nullptr));
        }

        std::vector<double> two_ptm_scr;
        for (size_t i = 0; i < ptm_pair_vec.size(); i++) {
            change_ptr1->setPtmPtr(ptm_pair_vec[i].first);
            change_ptr2->setPtmPtr(ptm_pair_vec[i].second);
            two_ptm_scr.push_back(getScr(prsm, true, mng_ptr));
        }

        return *std::max_element(two_ptm_scr.begin(), two_ptm_scr.end());
    }
}

double getScr(PrsmPtr & prsm, bool known, const LocalMngPtr& mng_ptr) {

    IonTypePtr n_ion_type_ptr =
        prsm->getDeconvMsPtrVec()[0]->getHeaderPtr()->getActivationPtr()->getNIonTypePtr();

    std::vector<double> theo_double =
        prsm->getProteoformPtr()->getBpSpecPtr()->getBreakPointMasses(n_ion_type_ptr);

    // theo_double includes a zero in the head
    std::sort(theo_double.begin(), theo_double.end());

    SpectrumSetPtr spec_set_ptr =
        getSpectrumSet(prsm->getDeconvMsPtrVec(), mng_ptr->prsm_para_ptr_->getSpParaPtr(),
                       prsm->getAdjustedPrecMass() );

    ExtendMsPtrVec ms_ptr_vec = prsm->getRefineMsPtrVec();

    // no zero in spec_double
    std::vector<double> spec_double;
    std::vector<double> spec_torlerance;
    for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
        for (size_t j = 0; j < ms_ptr_vec[i]->getPeakPtrVec().size(); j++) {
            spec_double.push_back(ms_ptr_vec[i]->getPeakPtr(j)->getMonoMass());
            spec_torlerance.push_back(ms_ptr_vec[i]->getPeakPtr(j)->getOrigTolerance());
        }
    }

    double scr = 0.0;

    ChangePtr change_ptr1 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
    ChangePtr change_ptr2 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[1];

    int left_pos1 = change_ptr1->getLeftBpPos(), right_pos1 = change_ptr1
                    ->getRightBpPos();

    int left_pos2 = change_ptr2->getLeftBpPos(), right_pos2 = change_ptr2
                    ->getRightBpPos();

    PtmPtr ptm1 = change_ptr1->getPtmPtr();
    PtmPtr ptm2 = change_ptr2->getPtmPtr();

    if (!known) {
        massAdjust(prsm, mng_ptr);
    }

    double mass1 = change_ptr1->getMassShift();
    double mass2 = change_ptr2->getMassShift();

    // len = theo_double.size() - 1
    int len = prsm->getProteoformPtr()->getLen();
    int num_match = prsm->getMatchPeakNum();
    int two_ptm_table[3][len + 1][num_match + 1];
    memset(two_ptm_table, 0, sizeof(int) * 3 * (len + 1) * (num_match + 1));
    two_ptm_table[0][0][0] = 1;

    double theo_mod;
    int n_match;

    for (int i = 1; i <= len; i++) {
        theo_mod = theo_double[i];
        n_match = getNumMatch(theo_mod, spec_double, spec_torlerance);
        for (int j = 0; j < num_match + 1; j++)
            two_ptm_table[0][i][j] = two_ptm_table[0][i - 1][j - n_match];
    }

    for (int i = 1; i < len + 1; i++) {
        theo_mod = theo_double[i] + mass1;
        n_match = getNumMatch(theo_mod, spec_double, spec_torlerance);
        for (int j = 0; j < num_match + 1; j++) {
            two_ptm_table[1][i][j] = two_ptm_table[0][i - 1][j - n_match]
                                     + two_ptm_table[1][i - 1][j - n_match];
        }
    }

    for (int i = 1; i < len + 1; i++) {
        theo_mod = theo_double[i] + mass1 + mass2;
        n_match = getNumMatch(theo_mod, spec_double, spec_torlerance);
        for (int j = 0; j < num_match + 1; j++) {
            two_ptm_table[2][i][j] = two_ptm_table[1][i - 1][j - n_match]
                                     + two_ptm_table[2][i - 1][j - n_match];
        }
    }

    for (int j = 0; j <= num_match; j++) {
        scr += two_ptm_table[2][len][j] * std::pow(mng_ptr->p1_, len - j)
               * std::pow(mng_ptr->p2_, j);
    }

    change_ptr1->setMassShift(mass1);
    change_ptr1->setLeftBpPos(left_pos1);
    change_ptr1->setRightBpPos(right_pos1);

    change_ptr2->setMassShift(mass2);
    change_ptr2->setLeftBpPos(left_pos2);
    change_ptr2->setRightBpPos(right_pos2);

    return scr;
}

void massAdjust(PrsmPtr& prsm, const LocalMngPtr& mng_ptr) {
    double ppm = mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();

    ChangePtr change_ptr1 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
    ChangePtr change_ptr2 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[1];

    double mass1 = change_ptr1->getMassShift();
    double mass2 = change_ptr2->getMassShift();

    PtmPtr ptm1 = change_ptr1->getPtmPtr();
    PtmPtr ptm2 = change_ptr2->getPtmPtr();

    if (std::abs(mass1 + mass2 - ptm1->getMonoMass() - ptm2->getMonoMass())
            < prsm->getOriPrecMass() * ppm) {
        double err = mass1 + mass2 - ptm1->getMonoMass() - ptm2->getMonoMass();
        change_ptr1->setMassShift(ptm1->getMonoMass() + err / 2);
        change_ptr2->setMassShift(ptm2->getMonoMass() + err / 2);
    } else if (mass1 + mass2 > ptm1->getMonoMass() + ptm2->getMonoMass()) {
        double err = mass1 + mass2 - ptm1->getMonoMass() - ptm2->getMonoMass() - 1;
        change_ptr1->setMassShift(ptm1->getMonoMass() + err / 2 + 1);
        change_ptr2->setMassShift(ptm2->getMonoMass() + err / 2);
    } else {
        double err = mass1 + mass2 - ptm1->getMonoMass() - ptm2->getMonoMass() + 1;
        change_ptr1->setMassShift(ptm1->getMonoMass() + err / 2 - 1);
        change_ptr2->setMassShift(ptm2->getMonoMass() + err / 2);
    }
}

void termAdjust(PrsmPtr& prsm, bool& ptm1_known, bool &ptm2_known, double& t_two_ptm_scr, const LocalMngPtr& mng_ptr) {

    int left_sup, right_sup;
    double ppm = mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
    ChangePtr change_ptr1 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
    ChangePtr change_ptr2 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[1];

    int ori_start = prsm->getProteoformPtr()->getStartPos();
    int ori_end = prsm->getProteoformPtr()->getEndPos();

    double ori_mass1 = change_ptr1->getMassShift();
    double ori_mass2 = change_ptr2->getMassShift();
    double mass1 = ori_mass1, mass2 = ori_mass2;

    int tmp;
    getSupPeakNum(prsm, change_ptr1, mng_ptr->min_mass_, left_sup, tmp);
    getSupPeakNum(prsm, change_ptr2, mng_ptr->min_mass_, tmp, right_sup);

    LOG_DEBUG("left_sup " << left_sup << " right_sup " << right_sup);

    int n_trunc_min, n_trunc_max, c_trunc_min, c_trunc_max;
    getNtermTruncRange(prsm, n_trunc_min, n_trunc_max, mng_ptr->max_ptm_mass_);
    getCtermTruncRange(prsm, c_trunc_min, c_trunc_max, mng_ptr->max_ptm_mass_);

    std::vector<bool> ptm1_known_vec, ptm2_known_vec;
    std::vector<double> t_scr;
    std::vector<int> c_vec;
    std::vector<int> n_vec;
    std::string n_seq, c_seq;

    int n_t = 0, n_f = 0, n_t_f = 0;

    if (left_sup <= LEFT_SUP_LIMIT && right_sup <= RIGHT_SUP_LIMIT) {
        for (int i = n_trunc_min; i <= n_trunc_max; i++) {
            for (int j = c_trunc_min; j <= c_trunc_max; j++) {
                if (i < 0) {
                    n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                            .substr(ori_start + i, -i);
                    mass1 = ori_mass1 - AcidFactory::getPeptideMass(n_seq);
                } else {
                    n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                            .substr(ori_start, i);
                    mass1 = ori_mass1 + AcidFactory::getPeptideMass(n_seq);
                }

                if (j >= 0) {
                    c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                            .substr(ori_end + 1, j);
                    mass2 = ori_mass2 - AcidFactory::getPeptideMass(c_seq);
                } else {
                    c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                            .substr(ori_end + 1 + j, -j);
                    mass2 = ori_mass2 + AcidFactory::getPeptideMass(c_seq);
                }

                n_vec.push_back(i);
                c_vec.push_back(j);

                change_ptr1->setMassShift(mass1);
                change_ptr2->setMassShift(mass2);
                prsm->getProteoformPtr()->setStartPos(ori_start + i);
                prsm->getProteoformPtr()->setEndPos(ori_end + j);

                if (std::abs(mass1) <= 1 + prsm->getOriPrecMass() * ppm
                        && std::abs(mass2) <= 1 + prsm->getOriPrecMass() * ppm) {
                    return;
                }

                bool is_known1 = PtmFactory::isKnown(mass1,
                                                     prsm->getOriPrecMass() * ppm);
                bool is_known2 = PtmFactory::isKnown(mass2,
                                                     prsm->getOriPrecMass() * ppm);
                ptm1_known_vec.push_back(is_known1);
                ptm2_known_vec.push_back(is_known2);
                t_scr.push_back(getScr(prsm, mng_ptr));
            }
        }

        for (size_t i = 0; i < ptm1_known_vec.size(); i++) {
            if (ptm1_known_vec[i] && ptm2_known_vec[i]) {
                n_t++;
            } else if (ptm1_known_vec[i] || ptm2_known_vec[i]) {
                n_t_f++;
            } else {
                n_f++;
            }
        }

        for (size_t i = 0; i < ptm1_known_vec.size(); i++) {
            if (ptm1_known_vec[i] && ptm2_known_vec[i]) {
                t_scr[i] = t_scr[i] * mng_ptr->theta_ * mng_ptr->theta_ / n_t;
            } else if (ptm1_known_vec[i] || ptm2_known_vec[i]) {
                t_scr[i] = t_scr[i] * mng_ptr->theta_ * (1 - mng_ptr->theta_) / n_t_f;
            } else {
                t_scr[i] = t_scr[i] * (1 - mng_ptr->theta_) * (1 - mng_ptr->theta_) / n_f;
            }
        }

        int idx = std::distance(t_scr.begin(),
                                std::max_element(t_scr.begin(), t_scr.end()));
        ptm1_known = ptm1_known_vec[idx];
        ptm2_known = ptm2_known_vec[idx];
        t_two_ptm_scr = t_scr[idx];
        prsm->getProteoformPtr()->setStartPos(ori_start + n_vec[idx]);
        prsm->getProteoformPtr()->setEndPos(ori_end + c_vec[idx]);

        if (n_vec[idx] < 0) {
            n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_start + n_vec[idx], -n_vec[idx]);
            mass1 = ori_mass1 - AcidFactory::getPeptideMass(n_seq);
        } else {
            n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_start, n_vec[idx]);
            mass1 = ori_mass1 + AcidFactory::getPeptideMass(n_seq);
        }

        if (c_vec[idx] >= 0) {
            c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_end + 1, c_vec[idx]);
            mass2 = ori_mass2 - AcidFactory::getPeptideMass(c_seq);
        } else {
            c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_end + 1 + c_vec[idx], -c_vec[idx]);
            mass2 = ori_mass2 + AcidFactory::getPeptideMass(c_seq);
        }
        change_ptr1->setMassShift(mass1);
        change_ptr2->setMassShift(mass2);

    } else if (left_sup <= LEFT_SUP_LIMIT) {
        for (int i = n_trunc_min; i <= n_trunc_max; i++) {
            n_vec.push_back(i);
            if (i < 0) {
                n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                        .substr(ori_start + i, -i);
                mass1 = ori_mass1 - AcidFactory::getPeptideMass(n_seq);
            } else {
                n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                        .substr(ori_start, i);
                mass1 = ori_mass1 + AcidFactory::getPeptideMass(n_seq);
            }

            change_ptr1->setMassShift(mass1);
            prsm->getProteoformPtr()->setStartPos(ori_start + i);

            if (std::abs(mass1) <= 1 + prsm->getOriPrecMass() * ppm) {
                return;
            }

            ptm1_known_vec.push_back(
                PtmFactory::isKnown(mass1, prsm->getOriPrecMass() * ppm));
            t_scr.push_back(getScr(prsm, mng_ptr));
        }

        for (size_t i = 0; i < ptm1_known_vec.size(); i++) {
            if (ptm1_known_vec[i] && ptm2_known) {
                n_t++;
            } else if (ptm1_known_vec[i] || ptm2_known) {
                n_t_f++;
            } else {
                n_f++;
            }
        }

        for (size_t i = 0; i < t_scr.size(); i++) {
            if (ptm1_known_vec[i] && ptm2_known) {
                t_scr[i] = t_scr[i] * mng_ptr->theta_ * mng_ptr->theta_ / n_t;
            } else if (ptm1_known_vec[i] || ptm2_known) {
                t_scr[i] = t_scr[i] * mng_ptr->theta_ * (1 - mng_ptr->theta_) / n_t_f;
            } else {
                t_scr[i] = t_scr[i] * (1 - mng_ptr->theta_) * (1 - mng_ptr->theta_) / n_f;
            }
        }
        int idx = std::distance(t_scr.begin(),
                                std::max_element(t_scr.begin(), t_scr.end()));

        ptm1_known = ptm1_known_vec[idx];
        t_two_ptm_scr = t_scr[idx];
        prsm->getProteoformPtr()->setStartPos(ori_start + n_vec[idx]);

        if (n_vec[idx] < 0) {
            n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_start + n_vec[idx], -n_vec[idx]);
            mass1 = ori_mass1 - AcidFactory::getPeptideMass(n_seq);
        } else {
            n_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_start, n_vec[idx]);
            mass1 = ori_mass1 + AcidFactory::getPeptideMass(n_seq);
        }
        change_ptr1->setMassShift(mass1);

    } else if (right_sup <= RIGHT_SUP_LIMIT) {
        LOG_DEBUG(c_trunc_min << " " << c_trunc_max);
        for (int i = c_trunc_min; i <= c_trunc_max; i++) {
            c_vec.push_back(i);

            if (i >= 0) {
                c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                        .substr(ori_end + 1, i);
                mass2 = ori_mass2 - AcidFactory::getPeptideMass(c_seq);
            } else {
                c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString()
                        .substr(ori_end + 1 + i, -i);
                mass2 = ori_mass2 + AcidFactory::getPeptideMass(c_seq);
            }

            change_ptr2->setMassShift(mass2);
            prsm->getProteoformPtr()->setEndPos(ori_end + i);

            if (std::abs(mass2) <= 1 + prsm->getOriPrecMass() * ppm) {
                return;
            }

            ptm2_known_vec.push_back(
                PtmFactory::isKnown(mass2, prsm->getOriPrecMass() * ppm));
            t_scr.push_back(getScr(prsm, mng_ptr));
        }

        for (size_t i = 0; i < ptm2_known_vec.size(); i++) {
            if (ptm1_known && ptm2_known_vec[i]) {
                n_t++;
            } else if (ptm1_known || ptm2_known_vec[i]) {
                n_t_f++;
            } else {
                n_f++;
            }
        }

        for (size_t i = 0; i < t_scr.size(); i++) {

            if (ptm1_known && ptm2_known_vec[i]) {
                t_scr[i] = t_scr[i] * mng_ptr->theta_ * mng_ptr->theta_ / n_t;
            } else if (ptm1_known || ptm2_known_vec[i]) {
                t_scr[i] = t_scr[i] * mng_ptr->theta_ * (1 - mng_ptr->theta_) / n_t_f;
            } else {
                t_scr[i] = t_scr[i] * (1 - mng_ptr->theta_) * (1 - mng_ptr->theta_) / n_f++;
            }
        }
        int idx = std::distance(t_scr.begin(),
                                std::max_element(t_scr.begin(), t_scr.end()));
        ptm2_known = ptm2_known_vec[idx];
        t_two_ptm_scr = t_scr[idx];
        prsm->getProteoformPtr()->setEndPos(ori_end + c_vec[idx]);

        if (c_vec[idx] >= 0) {
            c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_end + 1, c_vec[idx]);
            mass2 = ori_mass2 - AcidFactory::getPeptideMass(c_seq);
        } else {
            c_seq = prsm->getProteoformPtr()->getDbResSeqPtr()->toAcidString().substr(
                        ori_end + 1 + c_vec[idx], -c_vec[idx]);
            mass2 = ori_mass2 + AcidFactory::getPeptideMass(c_seq);
        }
        change_ptr2->setMassShift(mass2);
    } else {
        t_two_ptm_scr = getScr(prsm, mng_ptr);
        if (ptm1_known)
            t_two_ptm_scr *= mng_ptr->theta_;
        else
            t_two_ptm_scr *= (1 - mng_ptr->theta_);

        if (ptm2_known)
            t_two_ptm_scr *= mng_ptr->theta_;
        else
            t_two_ptm_scr *= (1 - mng_ptr->theta_);
    }
}

int getSplit(const PrsmPtr & prsm, const PtmPtr & ptm1, const PtmPtr & ptm2, const LocalMngPtr& mng_ptr) {
    IonTypePtr n_ion_type_ptr =
        prsm->getDeconvMsPtrVec()[0]->getHeaderPtr()->getActivationPtr()->getNIonTypePtr();

    std::vector<double> theo_double =
        prsm->getProteoformPtr()->getBpSpecPtr()->getBreakPointMasses(n_ion_type_ptr);

    std::sort(theo_double.begin(), theo_double.end());

    SpectrumSetPtr spec_set_ptr =
        getSpectrumSet(prsm->getDeconvMsPtrVec(),
                       mng_ptr->prsm_para_ptr_->getSpParaPtr(),
                       prsm->getAdjustedPrecMass());

    ExtendMsPtrVec ms_ptr_vec = prsm->getRefineMsPtrVec();

    std::vector<double> spec_double;//(ms_three_ptr->getPeakPtrVec().size());
    std::vector<double> spec_torlerance;//(spec_double.size());
    for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
        for (size_t j = 0; j < ms_ptr_vec[i]->getPeakPtrVec().size(); j++) {
            spec_double.push_back(ms_ptr_vec[i]->getPeakPtr(j)->getMonoMass());
            spec_torlerance.push_back(ms_ptr_vec[i]->getPeakPtr(j)->getMonoMass());
        }
    }

    ChangePtr change_ptr1 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
    ChangePtr change_ptr2 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[1];

    double mass1 = change_ptr1->getMassShift();
    double mass2 = change_ptr2->getMassShift();

    int len = prsm->getProteoformPtr()->getLen();
    int num_match = prsm->getMatchPeakNum();
    std::vector<double> split_scr;

    int two_ptm_table0[len + 1][num_match + 1];
    int two_ptm_table1[len + 1][num_match + 1];
    int two_ptm_table2[len + 1][num_match + 1];
    memset(two_ptm_table0, 0, sizeof(int) * (len + 1) * (num_match + 1));

    two_ptm_table0[0][0] = 1;

    double theo_mod;
    int n_match;
    for (int i = 1; i < len + 1; i++) {
        theo_mod = theo_double[i];
        n_match = getNumMatch(theo_mod, spec_double, spec_torlerance);
        for (int j = 0; j < num_match + 1; j++) {
            two_ptm_table0[i][j] = two_ptm_table0[i - 1][j - n_match];
        }
    }

    for (int i = 1; i < len; i++) {

        memset(two_ptm_table1, 0, sizeof(int) * (len + 1) * (num_match + 1));
        memset(two_ptm_table2, 0, sizeof(int) * (len + 1) * (num_match + 1));

        for (int j = 1; j <= i; j++) {
            theo_mod = theo_double[j] + mass1;
            n_match = getNumMatch(theo_mod, spec_double, spec_torlerance);
            for (int k = 1; k <= num_match; k++) {
                two_ptm_table1[j][k] = two_ptm_table0[j - 1][k - n_match]
                                       + two_ptm_table1[j - 1][k - n_match];
            }
        }

        for (int j = i + 1; j < len + 1; j++) {
            theo_mod = theo_double[j] + mass1 + mass2;
            n_match = getNumMatch(theo_mod, spec_double, spec_torlerance);
            for (int k = 1; k < num_match + 1; k++) {
                two_ptm_table2[j][k] = two_ptm_table1[j - 1][k - n_match]
                                       + two_ptm_table2[j - 1][k - n_match];
            }
        }

        double scr = 0.0;
        for (int j = 0; j <= num_match; j++) {
            scr += two_ptm_table2[len - 1][j] * std::pow(mng_ptr->p1_, len - j)
                   * std::pow(mng_ptr->p2_, j);
        }
        split_scr.push_back(scr);
    }

    return std::distance(split_scr.begin(), std::max_element(split_scr.begin(), split_scr.end()));
}

int getNumMatch(double theo_mod, const std::vector<double>& spec_double,
                const std::vector<double>& spec_torlerance) {

    int count = 0;
    for (size_t i = 0; i < spec_double.size(); i++) {
        if (std::abs(theo_mod - spec_double[i]) <= spec_torlerance[i]) {
            count++;
        }
    }
    return count;
}

} // namespace prot
