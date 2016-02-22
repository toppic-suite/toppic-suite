
#include "base/ptm.hpp"
#include "base/mod_util.hpp"
#include "base/algorithm.hpp"
#include "base/web_logger.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "spec/extend_ms_factory.hpp"
#include "prsm/prsm.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/peak_ion_pair_factory.hpp"
#include "spec/msalign_util.hpp"
#include "local_processor.hpp"
#include "local_util.hpp"
#include "local_anno.hpp"

namespace prot {

LocalProcessor::LocalProcessor(LocalMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  ppm_ = mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
  theta_ = mng_ptr->theta_;
  threshold_ = mng_ptr->threshold_;
  beta_ = mng_ptr->beta_;
  min_mass_ = mng_ptr->min_mass_;
  p1_ = mng_ptr->p1_;
  p2_ = mng_ptr->p2_;
  LocalUtil::init(mng_ptr_);
}

void LocalProcessor::process() {
  std::string spec_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = FileUtil::basename(spec_file_name) + "." + mng_ptr_->input_file_ext_;
  std::string output_file_name = FileUtil::basename(spec_file_name) + "." + mng_ptr_->output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();

  FastaIndexReaderPtr seq_reader(new FastaIndexReader(db_file_name));
  PrsmReader prsm_reader(input_file_name);
  ModPtrVec fix_mod_list = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_list);
  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  MsAlignReader sp_reader(spec_file_name, group_spec_num);
  SpectrumSetPtr spec_set_ptr;
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();

  int spectrum_num = MsAlignUtil::getSpNum (mng_ptr_->prsm_para_ptr_->getSpectrumFileName());
  int cnt = 0;
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    cnt += group_spec_num;
    if(spec_set_ptr->isValid()){
      int spec_id = spec_set_ptr->getSpecId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec 
            = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);

        if (prsm_ptr->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED) > 0)
          processOneSpectrum(prsm_ptr);

        writer.write(prsm_ptr);
        prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_list);
      }
    }
    std::cout << std::flush << "Localizaton is processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";
    WebLog::percentLog(cnt, spectrum_num, WebLog::LocalizationTime());
  }
  sp_reader.close();
  prsm_reader.close();
  writer.close();
  std::cout << std::endl;
}

void LocalProcessor::processOneSpectrum(PrsmPtr prsm) {
  if (prsm->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED) == 1) {
    double mass = prsm->getProteoformPtr()->getChangePtrVec(ChangeType::UNEXPECTED)[0]->getMassShift();
    // if the mass shift is between [-1, 1], we don't characterize it
    if (std::abs(mass) <= 1 + prsm->getAdjustedPrecMass() * ppm_) return;
    processOnePtm(prsm);  
  } else if (prsm->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED) == 2) {
    processTwoPtm(prsm);
  }
}

void LocalProcessor::processOnePtm(PrsmPtr prsm) {

  // the original number of matched fragment ions
  int ori_num_match_ion = LocalUtil::compNumPeakIonPairs(prsm->getProteoformPtr(), prsm->getRefineMsPtrVec());

  // we will get a nullptr if the mass shift can't be explained by a known variable ptm
  ProteoformPtr one_known_prsm = processOneKnown(prsm);

  // if it can be explained by one known ptm and the decrease in matched fragement
  // ions is small
  if (one_known_prsm != nullptr 
      && LocalUtil::compNumPeakIonPairs(one_known_prsm, prsm->getRefineMsPtrVec()) 
      > ori_num_match_ion - DESC_MATCH_LIMIT) {
    prsm->setProteoformPtr(one_known_prsm);
    return;
  }

  ProteoformPtr one_unknown_prsm = processOneUnknown(prsm);
  ProteoformPtr two_known_prsm = processTwoKnown(prsm);

  if (two_known_prsm != nullptr) {
    if (two_known_prsm->getChangePtrVec(ChangeType::UNEXPECTED)[0]->getLocalAnno()->getRawScr()
        > one_unknown_prsm->getChangePtrVec(ChangeType::UNEXPECTED)[0]->getLocalAnno()->getRawScr()) {
      if (LocalUtil::compNumPeakIonPairs(two_known_prsm, prsm->getRefineMsPtrVec()) 
          > ori_num_match_ion - DESC_MATCH_LIMIT) {
        prsm->setProteoformPtr(two_known_prsm);
        return;
      }
    }
  }

}

void LocalProcessor::processTwoPtm(PrsmPtr prsm) {
  int ori_num_match_ion = LocalUtil::compNumPeakIonPairs(prsm->getProteoformPtr(), prsm->getRefineMsPtrVec());

  ProteoformPtr two_known_prsm = processTwoKnown(prsm);

  if (two_known_prsm != nullptr 
      && LocalUtil::compNumPeakIonPairs(two_known_prsm, prsm->getRefineMsPtrVec()) 
      > ori_num_match_ion - DESC_MATCH_LIMIT) {
    prsm->setProteoformPtr(two_known_prsm);
    return;
  }

  ProteoformPtr two_unknown_prsm = processTwoUnknown(prsm);
  ProteoformPtr one_known_prsm = processOneKnown(prsm);

  if (one_known_prsm != nullptr) {
    // no ptm
    if (one_known_prsm->getChangePtrVec(ChangeType::UNEXPECTED)[0]->getLocalAnno() == nullptr) {
      if (LocalUtil::compNumPeakIonPairs(one_known_prsm, prsm->getRefineMsPtrVec()) 
          > ori_num_match_ion - DESC_MATCH_LIMIT) {
        prsm->setProteoformPtr(one_known_prsm);
        return;
      }
    } else if (two_unknown_prsm != nullptr) {
      double one_known_scr = 
          one_known_prsm->getChangePtrVec(ChangeType::UNEXPECTED)[0]->getLocalAnno()->getRawScr();
      double two_unknown_scr = 
          two_unknown_prsm->getChangePtrVec(ChangeType::UNEXPECTED)[0]->getLocalAnno()->getRawScr();
      if (one_known_scr > two_unknown_scr) {
        if (LocalUtil::compNumPeakIonPairs(one_known_prsm, prsm->getRefineMsPtrVec()) 
            > ori_num_match_ion - DESC_MATCH_LIMIT) {
          prsm->setProteoformPtr(one_known_prsm);
          return;
        }
      } else {
        if (LocalUtil::compNumPeakIonPairs(two_unknown_prsm, prsm->getRefineMsPtrVec()) 
            > ori_num_match_ion - DESC_MATCH_LIMIT)
          prsm->setProteoformPtr(two_unknown_prsm);
      }
    }
  }
}

// we will get a nullptr if the mass shift can't be explained by a variable ptm
ProteoformPtr LocalProcessor::processOneKnown(const PrsmPtr & prsm) {
  ChangePtrVec change_vec = 
      prsm->getProteoformPtr()->getChangePtrVec(ChangeType::UNEXPECTED);

  double mass = 0.0;
  for (size_t i = 0; i < change_vec.size(); i++) {
    mass += change_vec[i]->getMassShift();
  }

  PtmPtrVec ptm_vec = 
      LocalUtil::getPtmPtrVecByMass(mass, prsm->getAdjustedPrecMass() * ppm_);

  ProteoformPtr one_known_proteoform;

  ChangePtrVec new_change_vec = getExpectedChangeVec(prsm->getProteoformPtr());
  new_change_vec.push_back(geneUnexpectedChange(change_vec[0], mass));
  std::sort(new_change_vec.begin(), new_change_vec.end(), Change::cmpPosInc);
  one_known_proteoform = 
      std::make_shared<Proteoform>(prsm->getProteoformPtr()->getFastaSeqPtr(), 
                                   prsm->getProteoformPtr()->getProtModPtr(), 
                                   prsm->getProteoformPtr()->getStartPos(),
                                   prsm->getProteoformPtr()->getEndPos(), 
                                   prsm->getProteoformPtr()->getResSeqPtr(),
                                   new_change_vec);
  if (ptm_vec.size() == 0) {
    LocalUtil::onePtmTermAdjust(one_known_proteoform, 
                                prsm->getRefineMsPtrVec(),
                                mass, prsm->getAdjustedPrecMass() * ppm_);                                     
    if (std::abs(mass) < 1 + prsm->getAdjustedPrecMass() * ppm_){
      ChangePtr change_ptr = one_known_proteoform->getChangePtrVec(ChangeType::UNEXPECTED)[0];
      change_ptr->setLeftBpPos(std::max(change_vec[0]->getLeftBpPos()
                                        + prsm->getProteoformPtr()->getStartPos() 
                                        - one_known_proteoform->getStartPos(), 0));
      change_ptr->setRightBpPos(std::min(change_vec[0]->getRightBpPos() + prsm->getProteoformPtr()->getStartPos(),
                                         one_known_proteoform->getEndPos()) - one_known_proteoform->getStartPos());
      return one_known_proteoform;
    }
    ptm_vec = LocalUtil::getPtmPtrVecByMass(mass, prsm->getAdjustedPrecMass() * ppm_);
  }
  // even after adjusting N/C-termimals, still no explanation. return nullptr
  if (ptm_vec.size() == 0) {
    new_change_vec = getExpectedChangeVec(prsm->getProteoformPtr());
    for (size_t i = 0; i < new_change_vec.size(); i++) {
      int left = one_known_proteoform->getStartPos() + new_change_vec[i]->getLeftBpPos() - prsm->getProteoformPtr()->getStartPos();
      int right = one_known_proteoform->getStartPos() + new_change_vec[i]->getRightBpPos() - prsm->getProteoformPtr()->getStartPos();
      new_change_vec[i]->setLeftBpPos(left);
      new_change_vec[i]->setRightBpPos(right);
    }
    return nullptr;
  }

  double raw_scr;
  std::vector<double> scr_vec;
  ExtendMsPtrVec extend_ms_ptr_vec = prsm->getRefineMsPtrVec();

  LocalUtil::compOnePtmScr(one_known_proteoform, extend_ms_ptr_vec,
                           scr_vec, raw_scr, ptm_vec);
  // can be explained by a variable ptm, but no modifiable site
  if (scr_vec.size() == 0) return nullptr;

  int bgn, end;
  double conf;
  LocalUtil::scr_filter(scr_vec, bgn, end, conf, threshold_);
  if (bgn == -1) return nullptr;
  // it is known ptm, raw_scr * theta_; otherwise raw_scr * (1 - theta_)
  LocalAnnoPtr anno = std::make_shared<LocalAnno>(bgn, end, conf, scr_vec, raw_scr * theta_, ptm_vec[0]);
  one_known_proteoform->getChangePtrVec(ChangeType::UNEXPECTED)[0]->setLocalAnno(anno);
  return one_known_proteoform;
}

// no nullptr will be returned
ProteoformPtr LocalProcessor::processOneUnknown(const PrsmPtr & prsm) {
  ChangePtrVec change_vec = 
      prsm->getProteoformPtr()->getChangePtrVec(ChangeType::UNEXPECTED);

  double mass = 0.0;
  for (size_t i = 0; i < change_vec.size(); i++) {
    mass += change_vec[i]->getMassShift();
  }

  ProteoformPtr one_unknown_proteoform;

  ChangePtrVec new_change_vec = getExpectedChangeVec(prsm->getProteoformPtr());
  new_change_vec.push_back(geneUnexpectedChange(change_vec[0], mass));
  std::sort(new_change_vec.begin(), new_change_vec.end(), Change::cmpPosInc);
  one_unknown_proteoform = 
      std::make_shared<Proteoform>(prsm->getProteoformPtr()->getFastaSeqPtr(), 
                                   prsm->getProteoformPtr()->getProtModPtr(), 
                                   prsm->getProteoformPtr()->getStartPos(),
                                   prsm->getProteoformPtr()->getEndPos(), 
                                   prsm->getProteoformPtr()->getResSeqPtr(),
                                   new_change_vec);    

  LocalUtil::onePtmTermAdjust(one_unknown_proteoform, 
                              prsm->getRefineMsPtrVec(), mass,
                              prsm->getAdjustedPrecMass() * ppm_);

  double raw_scr;
  std::vector<double> scr_vec;
  PtmPtrVec ptm_vec;
  ExtendMsPtrVec extend_ms_ptr_vec = prsm->getRefineMsPtrVec();

  LocalUtil::compOnePtmScr(one_unknown_proteoform, extend_ms_ptr_vec,
                           scr_vec, raw_scr, ptm_vec);

  int bgn, end;
  double conf;
  LocalUtil::scr_filter(scr_vec, bgn, end, conf, threshold_);
  // it is known ptm, raw_scr * theta_; otherwise raw_scr * (1 - theta_)
  LocalAnnoPtr anno = std::make_shared<LocalAnno>(bgn, end, conf, scr_vec, raw_scr * (1 - theta_), nullptr);
  one_unknown_proteoform->getChangePtrVec(ChangeType::UNEXPECTED)[0]->setLocalAnno(anno);
  return one_unknown_proteoform;
}

// similar to processOneKnown, we might get a nullptr from this function
ProteoformPtr LocalProcessor::processTwoKnown(const PrsmPtr & prsm) {
  ChangePtrVec change_vec = prsm->getProteoformPtr()->getChangePtrVec(ChangeType::UNEXPECTED);

  ChangePtr change_ptr1 = geneUnexpectedChange(change_vec[0], change_vec[0]->getMassShift());
  ChangePtr change_ptr2;

  if (change_vec.size() == 1) {
    change_ptr2 = geneUnexpectedChange(change_vec[0], 0.0);
  } else {
    change_ptr2 = geneUnexpectedChange(change_vec[1], change_vec[1]->getMassShift());
  }

  double mass1 = change_ptr1->getMassShift();
  double mass2 = change_ptr2->getMassShift();

  PtmPairVec ptm_pair_vec = 
      LocalUtil::getPtmPairVecByMass(mass1, mass2, prsm->getAdjustedPrecMass() * ppm_);

  ProteoformPtr two_known_proteoform;

  ChangePtrVec new_change_vec = getExpectedChangeVec(prsm->getProteoformPtr());
  new_change_vec.push_back(change_ptr1);
  new_change_vec.push_back(change_ptr2);
  std::sort(new_change_vec.begin(), new_change_vec.end(), Change::cmpPosInc);
  two_known_proteoform = 
      std::make_shared<Proteoform>(prsm->getProteoformPtr()->getFastaSeqPtr(), 
                                   prsm->getProteoformPtr()->getProtModPtr(), 
                                   prsm->getProteoformPtr()->getStartPos(),
                                   prsm->getProteoformPtr()->getEndPos(), 
                                   prsm->getProteoformPtr()->getResSeqPtr(),
                                   new_change_vec);

  if (ptm_pair_vec.size() == 0) {
    LocalUtil::twoPtmTermAdjust(two_known_proteoform, prsm->getMatchPeakNum(),
                                prsm->getRefineMsPtrVec(), prsm->getAdjustedPrecMass(),
                                mass1, mass2);
    ptm_pair_vec = LocalUtil::getPtmPairVecByMass(mass1, mass2, prsm->getAdjustedPrecMass() * ppm_);
  }

  if (ptm_pair_vec.size() == 0) {
    new_change_vec = getExpectedChangeVec(prsm->getProteoformPtr());
    for (size_t i = 0; i < new_change_vec.size(); i++) {
      int left = two_known_proteoform->getStartPos() + new_change_vec[i]->getLeftBpPos() - prsm->getProteoformPtr()->getStartPos();
      int right = two_known_proteoform->getStartPos() + new_change_vec[i]->getRightBpPos() - prsm->getProteoformPtr()->getStartPos();
      new_change_vec[i]->setLeftBpPos(left);
      new_change_vec[i]->setRightBpPos(right);
    }
    return nullptr;
  }

  double raw_scr;
  ExtendMsPtrVec extend_ms_ptr_vec = prsm->getRefineMsPtrVec();

  LocalUtil::compTwoPtmScr(two_known_proteoform, prsm->getMatchPeakNum(), extend_ms_ptr_vec,
                           prsm->getAdjustedPrecMass(),
                           raw_scr, ptm_pair_vec);
  // can be explained by a variable ptm, but no modifiable site
  if (raw_scr == 0) return nullptr;

  PtmPtr ptm1 = ptm_pair_vec[0].first;
  PtmPtr ptm2 = ptm_pair_vec[0].second;
  raw_scr = raw_scr * theta_ * theta_ * (1 - beta_);
  std::vector<double> empty_scr_vec;
  LocalAnnoPtr anno1 = std::make_shared<LocalAnno>(0, 0, 0, empty_scr_vec, raw_scr, ptm1);
  two_known_proteoform->getChangePtrVec(ChangeType::UNEXPECTED)[0]->setLocalAnno(anno1);
  LocalAnnoPtr anno2 = std::make_shared<LocalAnno>(0, 0, 0, empty_scr_vec, raw_scr, ptm2);
  two_known_proteoform->getChangePtrVec(ChangeType::UNEXPECTED)[1]->setLocalAnno(anno2);
  LocalUtil::compSplitPoint(two_known_proteoform,prsm->getMatchPeakNum(),
                            extend_ms_ptr_vec, prsm->getAdjustedPrecMass());
  return two_known_proteoform;
}

ProteoformPtr LocalProcessor::processTwoUnknown(const PrsmPtr & prsm) {
  ChangePtrVec change_vec = prsm->getProteoformPtr()->getChangePtrVec(ChangeType::UNEXPECTED);

  ChangePtr change_ptr1 = geneUnexpectedChange(change_vec[0], change_vec[0]->getMassShift());
  ChangePtr change_ptr2;

  if (change_vec.size() == 1) {
    change_ptr2 = geneUnexpectedChange(change_vec[0], 0.0);
  } else {
    change_ptr2 = geneUnexpectedChange(change_vec[1], change_vec[1]->getMassShift());
  }

  double mass1 = change_ptr1->getMassShift();
  double mass2 = change_ptr2->getMassShift();
  PtmPtrVec ptm_vec1 = LocalUtil::getPtmPtrVecByMass(mass1, prsm->getAdjustedPrecMass() * ppm_);
  if (ptm_vec1.size() == 0) ptm_vec1.push_back(nullptr);

  PtmPtrVec ptm_vec2 = LocalUtil::getPtmPtrVecByMass(mass2, prsm->getAdjustedPrecMass() * ppm_);
  if (ptm_vec2.size() == 0) ptm_vec2.push_back(nullptr);

  PtmPairVec ptm_pair_vec;

  for (size_t i = 0; i < ptm_vec1.size(); i++) {
    for (size_t j = 0; j < ptm_vec2.size(); j++) {
      ptm_pair_vec.push_back(std::make_pair(ptm_vec1[i], ptm_vec2[j]));
    }
  }

  ChangePtrVec new_change_vec = getExpectedChangeVec(prsm->getProteoformPtr());
  new_change_vec.push_back(change_ptr1);
  new_change_vec.push_back(change_ptr2);
  std::sort(new_change_vec.begin(), new_change_vec.end(), Change::cmpPosInc);
  ProteoformPtr two_unknown_proteoform = 
      std::make_shared<Proteoform>(prsm->getProteoformPtr()->getFastaSeqPtr(), 
                                   prsm->getProteoformPtr()->getProtModPtr(), 
                                   prsm->getProteoformPtr()->getStartPos(),
                                   prsm->getProteoformPtr()->getEndPos(), 
                                   prsm->getProteoformPtr()->getResSeqPtr(),
                                   new_change_vec);

  double raw_scr;
  ExtendMsPtrVec extend_ms_ptr_vec = prsm->getRefineMsPtrVec();

  LocalUtil::compTwoPtmScr(two_unknown_proteoform, prsm->getMatchPeakNum(),
                           extend_ms_ptr_vec, prsm->getAdjustedPrecMass(),
                           raw_scr, ptm_pair_vec);

  PtmPtr ptm1 = ptm_pair_vec[0].first;
  if (ptm1 == nullptr)
    raw_scr = raw_scr * (1 - theta_);
  else
    raw_scr = raw_scr * theta_;

  PtmPtr ptm2 = ptm_pair_vec[0].second;
  if (ptm2 == nullptr)
    raw_scr = raw_scr * (1 - theta_);
  else
    raw_scr = raw_scr * theta_;

  std::vector<double> empty_scr_vec;
  LocalAnnoPtr anno1 = std::make_shared<LocalAnno>(0, 0, 0, empty_scr_vec, raw_scr, ptm1);
  two_unknown_proteoform->getChangePtrVec(ChangeType::UNEXPECTED)[0]->setLocalAnno(anno1);
  LocalAnnoPtr anno2 = std::make_shared<LocalAnno>(0, 0, 0, empty_scr_vec, raw_scr, ptm2);
  two_unknown_proteoform->getChangePtrVec(ChangeType::UNEXPECTED)[1]->setLocalAnno(anno2);

  LocalUtil::compSplitPoint(two_unknown_proteoform, prsm->getMatchPeakNum(),
                            extend_ms_ptr_vec, prsm->getAdjustedPrecMass());
  return two_unknown_proteoform;
}

}

