#include <string.h>
#include <boost/range/algorithm/remove_if.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "base/ptm.hpp"
#include "base/algorithm.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "prsm/peak_ion_pair.hpp"

#include "local_processor.hpp"
#include "local_util.hpp"

namespace prot {

LocalProcessor::LocalProcessor(LocalMngPtr& mng_ptr) {
    mng_ptr_ = mng_ptr;
    fai_ = NULL;
    ppm_ = mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
    p1_ = p2_ = 0.0;
    weight_ = mng_ptr->weight_;
    theta_ = mng_ptr->theta_;
    thread_ = mng_ptr->thread_;
    cysteine_protected_ = mng_ptr->cysteine_protected_;
    beta_ = mng_ptr->beta_;
    max_ptm_mass_ = mng_ptr->max_ptm_mass_;
    min_mass_ = mng_ptr->min_mass_;
}

void LocalProcessor::init() {
    //para_.resize(4);
    //std::fill(para_.begin(), para_.end(), 0.0);

    std::string spectrum_file_name =
        mng_ptr_->prsm_para_ptr_->getSpectrumFileName();

    std::string input_file_name = basename(spectrum_file_name) + "."
                                  + mng_ptr_->input_file_ext_;
    LOG_DEBUG("input file name " << input_file_name);

    ResiduePtrVec residue_ptr_vec = mng_ptr_->prsm_para_ptr_
                                    ->getFixModResiduePtrVec();

    prsm_ptrs_ = readAllPrsms(input_file_name,
                              mng_ptr_->prsm_para_ptr_->getSearchDbFileName(),
                              residue_ptr_vec);
    LOG_DEBUG("prsm loaded");
    addSpectrumPtrsToPrsms(prsm_ptrs_, mng_ptr_->prsm_para_ptr_);
    LOG_DEBUG("spectrum added");

    /*for (size_t i = 0; i < prsm_ptrs_.size(); i++) {*/
        //if (prsm_ptrs_[i]->getProteoformPtr()->getUnexpectedChangeNum() == 0) {

            //IonTypePtr n_ion_type_ptr =
                //prsm_ptrs_[i]->getDeconvMsPtrVec()[0]->getHeaderPtr()->getActivationPtr()
                //->getNIonTypePtr();

            //std::vector<double> theo_double =
                //prsm_ptrs_[i]->getProteoformPtr()->getBpSpecPtr()->getBreakPointMasses(
                    //n_ion_type_ptr);

            //std::vector<int> theo_int(theo_double.size());

            //std::transform(theo_double.begin(), theo_double.end(), theo_int.begin(),
            //[](double d) {
                //return std::round(d * SCALE_FACTOR);
            //});

            //std::sort(theo_int.begin(), theo_int.end());
            //auto last = std::unique(theo_int.begin(), theo_int.end());
            //theo_int.erase(last, theo_int.end());

            //std::vector<int> theo_int_ext;

            //for (size_t j = 0; j < theo_double.size(); j++) {
                //theo_int_ext.push_back(std::round((theo_double[j] + 42.01056) * SCALE_FACTOR));
                //theo_int_ext.push_back(std::round((theo_double[j] + 79.95682) * SCALE_FACTOR));
                //theo_int_ext.push_back(std::round((theo_double[j] + 14.01565) * SCALE_FACTOR));
                //theo_int_ext.push_back(std::round((theo_double[j] + 15.99491) * SCALE_FACTOR));
                //theo_int_ext.push_back(std::round((theo_double[j] + 114.042927) * SCALE_FACTOR));
                //theo_int_ext.push_back(std::round((theo_double[j] - 27.994915) * SCALE_FACTOR));
                //theo_int_ext.push_back(std::round((theo_double[j] - 18.0105) * SCALE_FACTOR));
            //}

            //// get unique
            //std::sort(theo_int.begin(), theo_int.end());
            //last = std::unique(theo_int.begin(), theo_int.end());
            //theo_int.erase(last, theo_int.end());

            //SpectrumSetPtr spec_set_ptr =
                //getSpectrumSet(prsm_ptrs_[i]->getDeconvMsPtrVec(),
                               //mng_ptr_->prsm_para_ptr_->getSpParaPtr(),
                               //prsm_ptrs_[i]->getAdjustedPrecMass());

            //PrmMsPtrVec prm = spec_set_ptr->getMsSixPtrVec();

            //std::vector<int> spec_int;
            //for (size_t j = 0; j < prm.size(); j++) {
                //for (size_t k = 0; k< prm[j]->getPeakPtrVec().size(); k++) {
                    //spec_int.push_back(std::round(prm[j]->getPeakPtr(k)->getMonoMass() * SCALE_FACTOR));
                //}
            //}

            //std::sort(spec_int.begin(), spec_int.end());
            //last = std::unique(spec_int.begin(), spec_int.end());
            //spec_int.erase(last, spec_int.end());

            //std::vector<int> common;

            //std::set_intersection(theo_int.begin(), theo_int.end(),
                                  //spec_int.begin(), spec_int.end(), std::back_inserter(common));

            //para_[0] += theo_int.size() - common.size();
            //para_[1] += theo_int_ext.size();
            //para_[2] += common.size();

            //common.clear();
            //std::set_intersection(theo_int_ext.begin(), theo_int_ext.end(),
                                  //spec_int.begin(), spec_int.end(), std::back_inserter(common));

            //para_[3] += common.size();
        //}
    //}

    //double p11, p01, p10, p00;

    //p11 = para_[2] / (para_[0] + para_[2]);
    //p01 = para_[0] / (para_[0] + para_[2]);
    //p10 = para_[3] * 7 / para_[1];
    //p00 = 1 - p10;

    //p1_ = p01 / p00;
    /*p2_ = p11 / p10;*/
    // const value of p1, p2 will be used in other datasets
    p1_ = 0.915258, p2_ = 21.1822;

    LOG_DEBUG("p1 = " << p1_ << " p2 = " << p2_ << " " << p2_ / p1_);
    mng_ptr_->p1_ = p1_;
    mng_ptr_->p2_ = p2_;
}

void LocalProcessor::process() {
    std::string spectrum_file_name =
        mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
    std::string output_file_name = basename(spectrum_file_name) + "."
                                   + mng_ptr_->output_file_ext_;
    PrsmWriter writer(output_file_name);
    int prsm_num = prsm_ptrs_.size();
    for (size_t i = 0; i < prsm_ptrs_.size(); i++) {
        processOneSpectrum(prsm_ptrs_[i], writer);
        std::cout << std::flush << "PTM localization is processing " << (i + 1)
                  << " of " << prsm_num << " prsms.\r";
    }
    writer.close();
}

void LocalProcessor::processOneSpectrum(PrsmPtr& prsm, PrsmWriter& writer) {

    if (prsm->getProteoformPtr()->getUnexpectedChangeNum() == 1) {
        processOnePtm(prsm);
        prsm->initMatchNum(mng_ptr_->min_mass_);
    } else if (prsm->getProteoformPtr()->getUnexpectedChangeNum() == 2) {
        processTwoPtm(prsm);
        prsm->initMatchNum(mng_ptr_->min_mass_);
    }
    writer.write(prsm);
}

void LocalProcessor::processOnePtm(PrsmPtr& prsm) {
    //LOG_DEBUG("One PTM, Scan " << prsm->getRefineMs()->getHeaderPtr()->getScansString());

    ChangePtr change_ptr =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];

    // original mass shift
    double ori_mass = change_ptr->getMassShift();

    int ori_start = prsm->getProteoformPtr()->getStartPos();
    int ori_end = prsm->getProteoformPtr()->getEndPos();

    if (std::abs(ori_mass) <= 1 + prsm->getOriPrecMass() * ppm_) {
        return;
    }

    bool ptm_known = PtmFactory::isKnown(ori_mass, prsm->getOriPrecMass() * ppm_);

    if (ptm_known) {
        processOneKnown(prsm);
        return;
    } else {  // truncation

        int left_sup = 0, right_sup = 0;
        getSupPeakNum(prsm, change_ptr, min_mass_, left_sup, right_sup);

        if (left_sup <= LEFT_SUP_LIMIT && prsm->getProteoformPtr()->isAcety()) {
            for (size_t j = 0; j < prsm->getProteoformPtr()->getChangePtrVec().size(); j++) {
                // the ptm might be the acetylation
                if (prsm->getProteoformPtr()->getChangePtrVec()[j]->getPtmPtr() != nullptr
                        && prsm->getProteoformPtr()->getChangePtrVec()[j]->
                        getPtmPtr()->getAbbrName() == "Acetylation") {
                    ori_mass += prsm->getProteoformPtr()->getChangePtrVec()[j]->getMassShift();
                    prsm->getProteoformPtr()->rmChangePtr(j);
                }
            }
        }

        change_ptr->setMassShift(ori_mass);

        if (std::abs(ori_mass) <= 1 + prsm->getOriPrecMass() * ppm_) {
            return;
        }

        double t_one_ptm_scr = 0.0;

        termAdjust(prsm, ptm_known, t_one_ptm_scr, mng_ptr_);

        if (prsm->getProteoformPtr()->getUnexpectedChangeNum() == 0)
            return;

        if (ptm_known) {
            processOneKnown(prsm);
            return;
        }

        PtmPairVec ptm_pair_vec = PtmFactory::getBasePtmPairByMass(
                                      ori_mass, 0, prsm->getOriPrecMass() * ppm_);

        if (ptm_pair_vec.size() > 0) {
            int t_start = prsm->getProteoformPtr()->getStartPos();
            int t_end = prsm->getProteoformPtr()->getEndPos();
            double t_mass = change_ptr->getMassShift();

            prsm->getProteoformPtr()->setStartPos(ori_start);
            prsm->getProteoformPtr()->setEndPos(ori_end);
            change_ptr->setMassShift(ori_mass);

            ChangePtr change_ptr2 = ChangePtr(
                                        new Change(0, 0, UNEXPECTED_CHANGE, 0, nullptr));

            prsm->getProteoformPtr()->addChangePtr(change_ptr2);

            std::vector<double> two_ptm_scr;
            for (size_t i = 0; i < ptm_pair_vec.size(); i++) {
                change_ptr->setPtmPtr(ptm_pair_vec[i].first);
                change_ptr2->setPtmPtr(ptm_pair_vec[i].second);
                two_ptm_scr.push_back(getScr(prsm, false, mng_ptr_));
            }

            int idx = std::distance(
                          two_ptm_scr.begin(),
                          std::max_element(two_ptm_scr.begin(), two_ptm_scr.end()));

            if (two_ptm_scr[idx] * (1 - beta_) * theta_ * theta_
                    > t_one_ptm_scr * beta_) {
                change_ptr->setPtmPtr(ptm_pair_vec[idx].first);
                change_ptr2->setPtmPtr(ptm_pair_vec[idx].second);

                massAdjust(prsm, mng_ptr_);
                processTwoKnown(prsm);

            } else {
                prsm->getProteoformPtr()->rmChangePtr(change_ptr2);
                change_ptr->setMassShift(t_mass);

                prsm->getProteoformPtr()->setStartPos(t_start);
                prsm->getProteoformPtr()->setEndPos(t_end);
                processOneUnknown(prsm);
            }
        } else {
            processOneUnknown(prsm);
        }

    }
}

void LocalProcessor::processOneKnown(PrsmPtr& prsm) {
    ChangePtr change_ptr =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];

    double mass_shift = change_ptr->getMassShift();

    std::vector<double> kScr, temp;
    std::vector<std::vector<double> > kScrVec;

    int n = prsm->getProteoformPtr()->getLen();

    // it is a known ptm, but we don't which one
    if (change_ptr->getPtmPtr() == nullptr) {
        PtmPtrVec ptmVec = PtmFactory::getBasePtmPtrByMass(
                               mass_shift, prsm->getOriPrecMass() * ppm_);
        for (size_t i = 0; i < ptmVec.size(); i++) {

            kScr.clear();
            for (int j = 0; j < n; j++) {
                if (modifiable(prsm->getProteoformPtr(), j, ptmVec[i],
                               cysteine_protected_)) {
                    change_ptr->setLeftBpPos(j);
                    change_ptr->setRightBpPos(j + 1);
                    int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                                   prsm->getRefineMsPtrVec(), min_mass_);
                    kScr.push_back(std::pow(p1_, n - match) * std::pow(p2_, match));
                } else {
                    kScr.push_back(0.0);
                }
            }
            kScrVec.push_back(kScr);
            temp.push_back(std::accumulate(kScr.begin(), kScr.end(), 0.0));
        }

        int idx = std::distance(temp.begin(),
                                std::max_element(temp.begin(), temp.end()));
        if (temp[idx] == 0) {
            processOneUnknown(prsm);
            return;
        }
        change_ptr->setPtmPtr(ptmVec[idx]);
        kScr = normalize(kScrVec[idx]);
    } else {
        for (int j = 0; j < n; j++) {
            if (modifiable(prsm->getProteoformPtr(), j, change_ptr->getPtmPtr(),
                           cysteine_protected_)) {
                change_ptr->setLeftBpPos(j);
                change_ptr->setRightBpPos(j + 1);
                int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                               prsm->getRefineMsPtrVec(), min_mass_);
                kScr.push_back(std::pow(p1_, n - match) * std::pow(p2_, match));
            } else {
                kScr.push_back(0.0);
            }
        }
        kScr = normalize(kScr);
    }

    int bgn, end;
    double conf;

    scr_filter(kScr, bgn, end, conf, thread_);
    change_ptr->setConf(conf);
    change_ptr->setScr(kScr);
    change_ptr->setLeftBpPos(bgn);
    change_ptr->setRightBpPos(end + 1);
    prsm->getProteoformPtr()->setSplit();
}

void LocalProcessor::processOneUnknown(PrsmPtr& prsm) {
    ChangePtr change_ptr =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];

    std::vector<double> ukScr;
    int n = prsm->getProteoformPtr()->getLen();

    for (int i = 0; i < n; i++) {
        if (modifiable(prsm->getProteoformPtr(), i, nullptr, cysteine_protected_)) {
            change_ptr->setLeftBpPos(i);
            change_ptr->setRightBpPos(i + 1);
            int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                           prsm->getRefineMsPtrVec(), min_mass_);
            ukScr.push_back(std::pow(p1_, n - match) * std::pow(p2_, match));
        } else {
            ukScr.push_back(0.0);
        }
    }

    ukScr = normalize(ukScr);

    int bgn, end;
    double conf;

    scr_filter(ukScr, bgn, end, conf, thread_);

    change_ptr->setPtmPtr(nullptr);
    change_ptr->setConf(conf);
    change_ptr->setScr(ukScr);
    change_ptr->setLeftBpPos(bgn);
    change_ptr->setRightBpPos(end + 1);
    prsm->getProteoformPtr()->setSplit();

}

void LocalProcessor::processTwoPtm(PrsmPtr& prsm) {

    //LOG_DEBUG("Two PTMs, Scan " << prsm->getRefineMs()->getHeaderPtr()->getScansString());

    ChangePtr change_ptr1 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
    ChangePtr change_ptr2 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[1];

    double mass1 = change_ptr1->getMassShift();
    double mass2 = change_ptr2->getMassShift();

    PtmPtr ptm1 = change_ptr1->getPtmPtr();
    PtmPtr ptm2 = change_ptr2->getPtmPtr();

    int ori_start = prsm->getProteoformPtr()->getStartPos();
    int ori_end = prsm->getProteoformPtr()->getEndPos();

    int left_sup, right_sup, tmp;

    getSupPeakNum(prsm, change_ptr1, min_mass_, left_sup, tmp);
    getSupPeakNum(prsm, change_ptr2, min_mass_, tmp, right_sup);

    bool ptm1_known = PtmFactory::isKnown(mass1, prsm->getOriPrecMass() * ppm_);
    bool ptm2_known = PtmFactory::isKnown(mass2, prsm->getOriPrecMass() * ppm_);

    if (ptm1_known && ptm2_known) {
        processTwoKnown(prsm);
        return;
    }

    double t_two_ptm_scr;
    termAdjust(prsm, ptm1_known, ptm2_known, t_two_ptm_scr, mng_ptr_);

    if (ptm1_known && ptm2_known) {
        processTwoKnown(prsm);
        return;
    }

    t_two_ptm_scr = t_two_ptm_scr * (1 - beta_);

    int t_start = prsm->getProteoformPtr()->getStartPos();
    int t_end = prsm->getProteoformPtr()->getEndPos();
    double t_mass1 = change_ptr1->getMassShift();
    double t_mass2 = change_ptr2->getMassShift();
    PtmPtr t_ptm1 = change_ptr1->getPtmPtr();
    PtmPtr t_ptm2 = change_ptr2->getPtmPtr();

    prsm->getProteoformPtr()->setStartPos(ori_start);
    prsm->getProteoformPtr()->setEndPos(ori_end);
    change_ptr1->setMassShift(mass1);
    change_ptr2->setMassShift(mass2);

    PtmPairVec ptm_pair_vec = PtmFactory::getBasePtmPairByMass(
                                  mass1, mass2, prsm->getOriPrecMass() * ppm_);

    int ptm_pair_vec_idx = 0;

    double two_ptm_scr = 0.0, one_ptm_scr = 0.0;

    if (ptm_pair_vec.size() > 0) {

        std::vector<double> two_ptm_scr_vec;

        for (size_t i = 0; i < ptm_pair_vec.size(); i++) {
            change_ptr1->setPtmPtr(ptm_pair_vec[i].first);
            change_ptr2->setPtmPtr(ptm_pair_vec[i].second);
            two_ptm_scr_vec.push_back(getScr(prsm, false, mng_ptr_));
        }

        two_ptm_scr = *std::max_element(two_ptm_scr_vec.begin(),
                                        two_ptm_scr_vec.end());

        two_ptm_scr = two_ptm_scr * theta_ * theta_ * (1 - beta_);

        LOG_DEBUG(two_ptm_scr);

        ptm_pair_vec_idx = std::distance(
                               two_ptm_scr_vec.begin(),
                               std::max_element(two_ptm_scr_vec.begin(), two_ptm_scr_vec.end()));
    }

    change_ptr1->setPtmPtr(nullptr);
    change_ptr2->setPtmPtr(nullptr);

    PtmPtrVec ptmVec = PtmFactory::getBasePtmPtrByMass(
                           mass1 + mass2, prsm->getOriPrecMass() * ppm_);

    prsm->getProteoformPtr()->rmChangePtr(change_ptr2);

    if (ptmVec.size() > 0) {
        change_ptr1->setMassShift(mass1 + mass2);
        std::vector<double> one_ptm_scr_vec;
        int n = prsm->getProteoformPtr()->getLen();
        for (size_t i = 0; i < ptmVec.size(); i++) {
            one_ptm_scr = 0.0;
            for (int j = 0; j < n; j++) {
                if (modifiable(prsm->getProteoformPtr(), j, ptmVec[i],
                               cysteine_protected_)) {
                    int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                                   prsm->getRefineMsPtrVec(), min_mass_);
                    one_ptm_scr += std::pow(p1_, n - match) * std::pow(p2_, match);
                }
            }
            one_ptm_scr_vec.push_back(one_ptm_scr);
        }
        one_ptm_scr = *std::max_element(one_ptm_scr_vec.begin(),
                                        one_ptm_scr_vec.end());
        one_ptm_scr = one_ptm_scr * beta_ * theta_;
    }

    if (t_two_ptm_scr == 0.0 && two_ptm_scr == 0.0 && one_ptm_scr == 0.0) {
        prsm->getProteoformPtr()->setStartPos(ori_start);
        prsm->getProteoformPtr()->setEndPos(ori_end);
        change_ptr1->setMassShift(mass1);
        change_ptr2->setMassShift(mass2);
        change_ptr1->setPtmPtr(ptm1);
        change_ptr2->setPtmPtr(ptm2);
        prsm->getProteoformPtr()->addChangePtr(change_ptr2);
        processTwoUnknown(prsm, ptm1_known, ptm2_known);
    }

    LOG_DEBUG(t_two_ptm_scr << " " << two_ptm_scr << " " << one_ptm_scr);

    if (one_ptm_scr > t_two_ptm_scr && one_ptm_scr > two_ptm_scr) {
        processOneKnown(prsm);
    } else if (two_ptm_scr > t_two_ptm_scr && two_ptm_scr > one_ptm_scr) {
        change_ptr1->setPtmPtr(ptm_pair_vec[ptm_pair_vec_idx].first);
        change_ptr2->setPtmPtr(ptm_pair_vec[ptm_pair_vec_idx].second);
        change_ptr1->setMassShift(t_mass1);
        change_ptr2->setMassShift(t_mass2);
        prsm->getProteoformPtr()->addChangePtr(change_ptr2);
        massAdjust(prsm, mng_ptr_);
        processTwoKnown(prsm);
    } else if (t_two_ptm_scr > two_ptm_scr && t_two_ptm_scr > one_ptm_scr) {
        prsm->getProteoformPtr()->setStartPos(t_start);
        prsm->getProteoformPtr()->setEndPos(t_end);
        change_ptr1->setMassShift(t_mass1);
        change_ptr2->setMassShift(t_mass2);
        change_ptr1->setPtmPtr(t_ptm1);
        change_ptr2->setPtmPtr(t_ptm2);
        prsm->getProteoformPtr()->addChangePtr(change_ptr2);
        processTwoUnknown(prsm, ptm1_known, ptm2_known);
    }
}

void LocalProcessor::processTwoKnown(PrsmPtr& prsm) {

    ChangePtr change_ptr1 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
    ChangePtr change_ptr2 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[1];

    int left_pos2 = change_ptr2->getLeftBpPos(), right_pos2 = change_ptr2
                    ->getRightBpPos();

    double mass1 = change_ptr1->getMassShift();
    double mass2 = change_ptr2->getMassShift();

    if (change_ptr1->getPtmPtr() == nullptr
            || change_ptr2->getPtmPtr() == nullptr) {
        PtmPtrVec ptmVec1 = PtmFactory::getBasePtmPtrByMass(
                                mass1, prsm->getOriPrecMass() * ppm_);
        PtmPtrVec ptmVec2 = PtmFactory::getBasePtmPtrByMass(
                                mass2, prsm->getOriPrecMass() * ppm_);

        PtmPairVec ptm_pair_vec;

        for (size_t i = 0; i < ptmVec1.size(); i++) {
            for (size_t j = 0; j < ptmVec2.size(); j++) {
                ptm_pair_vec.push_back(std::make_pair(ptmVec1[i], ptmVec2[j]));
            }
        }

        std::vector<double> two_ptm_scr;
        for (size_t i = 0; i < ptm_pair_vec.size(); i++) {
            change_ptr1->setPtmPtr(ptm_pair_vec[i].first);
            change_ptr2->setPtmPtr(ptm_pair_vec[i].second);
            two_ptm_scr.push_back(getScr(prsm, true, mng_ptr_));
        }

        int idx = std::distance(
                      two_ptm_scr.begin(),
                      std::max_element(two_ptm_scr.begin(), two_ptm_scr.end()));

        change_ptr1->setPtmPtr(ptm_pair_vec[idx].first);
        change_ptr2->setPtmPtr(ptm_pair_vec[idx].second);
    }

    PtmPtr ptm1 = change_ptr1->getPtmPtr();
    PtmPtr ptm2 = change_ptr2->getPtmPtr();

    int len = prsm->getProteoformPtr()->getLen();

    int split = getSplit(prsm, ptm1, ptm2, mng_ptr_);

    prsm->getProteoformPtr()->setSplit(split);

    std::vector<double> ptm1_scr;

    for (int i = 0; i <= split; i++) {
        change_ptr1->setLeftBpPos(i);
        change_ptr1->setRightBpPos(i + 1);

        if (modifiable(prsm->getProteoformPtr(), i, ptm1, cysteine_protected_)) {
            int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                           prsm->getRefineMsPtrVec(), min_mass_);
            ptm1_scr.push_back(std::pow(p1_, split - match) * std::pow(p2_, match));
        } else {
            ptm1_scr.push_back(0.0);
        }
    }

    ptm1_scr = normalize(ptm1_scr);

    int bgn, end;
    double conf;

    scr_filter(ptm1_scr, bgn, end, conf, thread_);

    change_ptr1->setConf(conf);
    change_ptr1->setScr(ptm1_scr);
    change_ptr1->setLeftBpPos(bgn);
    change_ptr1->setRightBpPos(end + 1);

    change_ptr2->setLeftBpPos(left_pos2);
    change_ptr2->setRightBpPos(right_pos2);

    std::vector<double> ptm2_scr;

    for (int i = split + 1; i < len; i++) {
        change_ptr2->setLeftBpPos(i);
        change_ptr2->setRightBpPos(i + 1);

        if (modifiable(prsm->getProteoformPtr(), i, ptm2, cysteine_protected_)) {
            int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                           prsm->getRefineMsPtrVec(), min_mass_);
            ptm2_scr.push_back(std::pow(p1_, len - match) * std::pow(p2_, match));
        } else {
            ptm2_scr.push_back(0.0);
        }
    }

    ptm2_scr = normalize(ptm2_scr);

    scr_filter(ptm2_scr, bgn, end, conf, thread_);

    change_ptr2->setConf(conf);
    change_ptr2->setScr(ptm2_scr);
    change_ptr2->setLeftBpPos(split + bgn + 1);
    change_ptr2->setRightBpPos(split + end + 2);
    change_ptr2->setSplit(split + 1);

}

void LocalProcessor::processTwoUnknown(PrsmPtr& prsm, bool ptm1_known, bool ptm2_known) {
    ChangePtr change_ptr1 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
    ChangePtr change_ptr2 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[1];
    PtmPtr ptm1 = change_ptr1->getPtmPtr();
    PtmPtr ptm2 = change_ptr2->getPtmPtr();

    int left_pos2 = change_ptr2->getLeftBpPos();
    int right_pos2 = change_ptr2->getRightBpPos();

    int len = prsm->getProteoformPtr()->getLen();

    int split = getSplit(prsm, ptm1, ptm2, mng_ptr_);

    prsm->getProteoformPtr()->setSplit(split);

    std::vector<double> ptm1_scr;
    std::vector<double> ukScr;

    for (int i = 0; i <= split; i++) {
        change_ptr1->setLeftBpPos(i);
        change_ptr1->setRightBpPos(i + 1);

        int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                       prsm->getRefineMsPtrVec(), min_mass_);
        if (modifiable(prsm->getProteoformPtr(), i, nullptr, cysteine_protected_)) {
            ukScr.push_back(std::pow(p1_, split - match) * std::pow(p2_, match));
        } else {
            ukScr.push_back(0.0);
        }
    }

    ukScr = normalize(ukScr);

    if (ptm1_known) {
        for (int i = 0; i <= split; i++) {
            change_ptr1->setLeftBpPos(i);
            change_ptr1->setRightBpPos(i + 1);

            if (modifiable(prsm->getProteoformPtr(), i, change_ptr1->getPtmPtr(),
                           cysteine_protected_)) {
                int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                               prsm->getRefineMsPtrVec(), min_mass_);
                ptm1_scr.push_back(std::pow(p1_, split - match) * std::pow(p2_, match));
            } else {
                ptm1_scr.push_back(0.0);
            }
        }
        ptm1_scr = normalize(ptm1_scr);
    } else {
        ptm1_scr = ukScr;
    }

    int bgn, end;
    double conf;

    scr_filter(ptm1_scr, bgn, end, conf, thread_);

    change_ptr1->setConf(conf);
    change_ptr1->setScr(ptm1_scr);
    change_ptr1->setLeftBpPos(bgn);
    change_ptr1->setRightBpPos(end + 1);

    change_ptr2->setLeftBpPos(left_pos2);
    change_ptr2->setRightBpPos(right_pos2);

    std::vector<double> ptm2_scr;
    ukScr.clear();

    for (int i = split + 1; i < len; i++) {
        change_ptr2->setLeftBpPos(i);
        change_ptr2->setRightBpPos(i + 1);

        if (modifiable(prsm->getProteoformPtr(), i, nullptr, cysteine_protected_)) {
            int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                           prsm->getRefineMsPtrVec(), min_mass_);
            ukScr.push_back(std::pow(p1_, len - match) * std::pow(p2_, match));
        } else {
            ukScr.push_back(0.0);
        }
    }

    ukScr = normalize(ukScr);

    if (ptm2_known) {
        for (int i = split + 1; i < len; i++) {
            change_ptr2->setLeftBpPos(i);
            change_ptr2->setRightBpPos(i + 1);

            if (modifiable(prsm->getProteoformPtr(), i, change_ptr2->getPtmPtr(),
                           cysteine_protected_)) {
                int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                               prsm->getRefineMsPtrVec(), min_mass_);
                ptm2_scr.push_back(std::pow(p1_, len - match) * std::pow(p2_, match));
            } else {
                ptm2_scr.push_back(0.0);
            }
        }
        ptm2_scr = normalize(ptm2_scr);
    } else {
        ptm2_scr = ukScr;
    }

    scr_filter(ptm2_scr, bgn, end, conf, thread_);

    change_ptr2->setConf(conf);
    change_ptr2->setScr(ptm2_scr);
    change_ptr2->setLeftBpPos(split + bgn + 1);
    change_ptr2->setRightBpPos(split + end + 2);
    change_ptr2->setSplit(split + 1);

}

}
