#include <fstream>
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
    for (size_t i = 0; i < prsm_ptrs_.size(); i++) {
        processOneSpectrum(prsm_ptrs_[i], writer);
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
    LOG_DEBUG("One PTM, Scan " << prsm->getRefineMs()->getHeaderPtr()->getScansString());

    ChangePtr change_ptr =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];

    double mass = change_ptr->getMassShift();

    int n = prsm->getProteoformPtr()->getLen();

    double one_ptm_scr = 0.0;

    int count = 0;

    PtmPtrVec ptmVec = PtmFactory::getBasePtmPtrByMass(
                           mass, prsm->getOriPrecMass() * ppm_);
    for (size_t i = 0; i < ptmVec.size(); i++) {

        if (ptmVec[i]->getName() != "Dimethylation")
            continue;

        for (int j = 0; j < n; j++) {
            count++;
            //if (modifiable(prsm->getProteoformPtr(), j, ptmVec[i],
            //               cysteine_protected_)) {
                change_ptr->setLeftBpPos(j);
                change_ptr->setRightBpPos(j + 1);
                int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                               prsm->getRefineMs(), min_mass_);
                one_ptm_scr += std::pow(p1_, n - match) * std::pow(p2_, match);
            //}
        }

    }

    one_ptm_scr = one_ptm_scr / count;

    PtmPairVec ptm_pair_vec = PtmFactory::getBasePtmPairByMass(mass, 0, prsm->getOriPrecMass() * ppm_);

    if (ptm_pair_vec.size() > 0) {

        ChangePtr change_ptr2 = ChangePtr(new Change(0, 0, UNEXPECTED_CHANGE, 0, nullptr));

        prsm->getProteoformPtr()->addChangePtr(change_ptr2);

        double two_ptm_scr = 0.0;

        count = 0;

        for (size_t i = 0; i < ptm_pair_vec.size(); i++) {
            if (ptm_pair_vec[i].first->getName() != "Methylation")
                continue;

            change_ptr->setPtmPtr(ptm_pair_vec[i].first);
            change_ptr2->setPtmPtr(ptm_pair_vec[i].second);
            for (int j = 0; j < n - 1; j++) {
                for (int k = j + 1; k < n; k++) {
                    count++;
                    //if (modifiable(prsm->getProteoformPtr(), j, ptm_pair_vec[i].first, cysteine_protected_)
                    //        && modifiable(prsm->getProteoformPtr(), k, ptm_pair_vec[i].second, cysteine_protected_) ) {

                        change_ptr->setLeftBpPos(j);
                        change_ptr->setRightBpPos(j + 1);
                        change_ptr2->setLeftBpPos(k);
                        change_ptr2->setRightBpPos(k + 1);

                        int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                                       prsm->getRefineMs(), min_mass_);

                        two_ptm_scr += std::pow(p1_, n - match) * std::pow(p2_, match);
                    //}
                }
            }
        }

        two_ptm_scr = two_ptm_scr / count;

        if(one_ptm_scr > 0 && two_ptm_scr > 0){
            std::ofstream outfile;
            outfile.open("one_two_score.txt", std::ios_base::app);
            outfile << one_ptm_scr << " " << two_ptm_scr << std::endl;
            outfile.close();
        }
    }

}

void LocalProcessor::processOneKnown(PrsmPtr& prsm) {}

void LocalProcessor::processOneUnknown(PrsmPtr& prsm) {}

void LocalProcessor::processTwoPtm(PrsmPtr& prsm) {

    LOG_DEBUG("Two PTMs, Scan " << prsm->getRefineMs()->getHeaderPtr()->getScansString());

    ChangePtr change_ptr1 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[0];
    ChangePtr change_ptr2 =
        prsm->getProteoformPtr()->getUnexpectedChangePtrVec()[1];

    double mass1 = change_ptr1->getMassShift();
    double mass2 = change_ptr2->getMassShift();

    std::cout << mass1 << " " << mass2 << std::endl;

    PtmPairVec ptm_pair_vec = PtmFactory::getBasePtmPairByMass(mass1, mass2, prsm->getOriPrecMass() * ppm_);

    double two_ptm_scr = 0.0, one_ptm_scr = 0.0;

    int n = prsm->getProteoformPtr()->getLen();

    int count = 0;

    if (ptm_pair_vec.size() > 0) {
        for (size_t i = 0 ; i < ptm_pair_vec.size(); i++) {

            if (ptm_pair_vec[i].first->getName() != "Methylation"
                || ptm_pair_vec[i].second->getName() != "Methylation")
                continue;

            change_ptr1->setPtmPtr(ptm_pair_vec[i].first);
            change_ptr2->setPtmPtr(ptm_pair_vec[i].second);
            count = 0;
            for (int j = 0; j < n - 1; j++) {
                for (int k = j + 1; k < n; k++) {

                    count++;

                    //if (modifiable(prsm->getProteoformPtr(), j, ptm_pair_vec[i].first, cysteine_protected_)
                    //        && modifiable(prsm->getProteoformPtr(), k, ptm_pair_vec[i].second, cysteine_protected_) ) {

                        change_ptr1->setLeftBpPos(j);
                        change_ptr1->setRightBpPos(j + 1);
                        change_ptr2->setLeftBpPos(k);
                        change_ptr2->setRightBpPos(k + 1);

                        int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                                       prsm->getRefineMs(), min_mass_);
                        two_ptm_scr += std::pow(p1_, n - match) * std::pow(p2_, match);
                    //}
                }
            }
        }
    }

    two_ptm_scr = two_ptm_scr / count;

    change_ptr1->setPtmPtr(nullptr);
    change_ptr2->setPtmPtr(nullptr);

    PtmPtrVec ptmVec = PtmFactory::getBasePtmPtrByMass(mass1 + mass2, prsm->getOriPrecMass() * ppm_);

    prsm->getProteoformPtr()->rmChangePtr(change_ptr2);

    if (ptmVec.size() > 0) {
        change_ptr1->setMassShift(mass1 + mass2);
        std::vector<double> one_ptm_scr_vec;
        int n = prsm->getProteoformPtr()->getLen();
        for (size_t i = 0; i < ptmVec.size(); i++) {

            if (ptmVec[i]->getName() != "Dimethylation")
                continue;

            count = 0;
            for (int j = 0; j < n; j++) {
                count++;
                //if (modifiable(prsm->getProteoformPtr(), j, ptmVec[i],cysteine_protected_)) {
                    change_ptr1->setLeftBpPos(j);
                    change_ptr1->setRightBpPos(j + 1);
                    int match = getNumPeakIonPairs(prsm->getProteoformPtr(),
                                                   prsm->getRefineMs(), min_mass_);
                    one_ptm_scr += std::pow(p1_, n - match) * std::pow(p2_, match);
                //}
            }
        }
    }

    one_ptm_scr = one_ptm_scr / count;
    if (one_ptm_scr > 0 && two_ptm_scr > 0){
        std::ofstream outfile;
        outfile.open("one_two_score.txt", std::ios_base::app);
        outfile << one_ptm_scr << " " << two_ptm_scr << std::endl;
        outfile.close();
    }
}

void LocalProcessor::processTwoKnown(PrsmPtr& prsm) {}

void LocalProcessor::processTwoUnknown(PrsmPtr& prsm, bool ptm1_known, bool ptm2_known) {}

}
