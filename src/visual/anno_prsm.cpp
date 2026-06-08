//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include "visual/anno_prsm.hpp"

#include "common/util/logger.hpp"
#include "ms/spec/peak_util.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "visual/anno_proteoform.hpp"

#include <iostream>

namespace toppic {

namespace anno_prsm {

void addPrsmHeader(XmlDOMDocument* xml_doc, XmlDOMElement element,
                   const PrsmPtr &prsm_ptr, const PrsmViewMngPtr &mng_ptr) {
  std::string str = std::to_string(prsm_ptr->getPrsmId());
  xml_doc->addElement(element, "prsm_id", str.c_str());
  if (prsm_ptr->getExpectedValuePtr().get() != nullptr) {
    str = str_util::evalueToString(prsm_ptr->getExpectedValuePtr()->getPValue(), 
                                   mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "p_value", str.c_str());
  } else {
    xml_doc->addElement(element, "p_value", "N/A");
  }

  if (prsm_ptr->getExpectedValuePtr().get() != nullptr) {
    str = str_util::evalueToString(prsm_ptr->getExpectedValuePtr()->getEValue(), 
                                   mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "e_value", str.c_str());
  } else {
    xml_doc->addElement(element, "e_value", "N/A");
  }

  double fdr = prsm_ptr->getFdr();

  if (fdr >= 0) {
    str = str_util::evalueToString(prsm_ptr->getFdr(), mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "fdr", str.c_str());
  } else {
    xml_doc->addElement(element, "fdr", "N/A");
  }

  str = std::to_string(static_cast<int>(prsm_ptr->getMatchFragNum()));

  xml_doc->addElement(element, "matched_fragment_number", str.c_str());

  str = std::to_string(static_cast<int>(prsm_ptr->getMatchPeakNum()));

  xml_doc->addElement(element, "matched_peak_number", str.c_str());
}

void addPrsmHeaderBrief(XmlDOMDocument* xml_doc, XmlDOMElement element,
                   const PrsmPtr &prsm_ptr, const PrsmViewMngPtr &mng_ptr) {
  std::string str = std::to_string(prsm_ptr->getPrsmId());
  xml_doc->addElement(element, "prsm_id", str.c_str());
}

void addMsHeader(XmlDOMDocument* xml_doc, XmlDOMElement ms_element, 
                 const PrsmPtr &prsm_ptr, const PrsmViewMngPtr &mng_ptr) {
  XmlDOMElement ms_header_element = xml_doc->addElement(ms_element, "ms_header");
  xml_doc->addElement(ms_header_element, "spectrum_file_name", prsm_ptr->getFileName().c_str());
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  std::string ms1_ids, ms2_ids;
  std::string ms1_scans, ms2_scans;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    MsHeaderPtr header_ptr = deconv_ms_ptr_vec[i]->getMsHeaderPtr();
    ms1_ids = ms1_ids + std::to_string(header_ptr->getMsOneId()) + " ";
    ms1_scans = ms1_scans + std::to_string(header_ptr->getMsOneScan()) + " ";
    ms2_ids = ms2_ids + std::to_string(header_ptr->getSpecId()) + " ";
    ms2_scans = ms2_scans + header_ptr->getScansString() + " ";
  }
  str_util::trim(ms1_ids);
  str_util::trim(ms1_scans);
  str_util::trim(ms2_ids);
  str_util::trim(ms2_scans);
  xml_doc->addElement(ms_header_element, "ms1_ids", ms1_ids.c_str());
  xml_doc->addElement(ms_header_element, "ms1_scans", ms1_scans.c_str());
  xml_doc->addElement(ms_header_element, "ids", ms2_ids.c_str());
  xml_doc->addElement(ms_header_element, "scans", ms2_scans.c_str());

  if (deconv_ms_ptr_vec.size() > 0) {
    int pos = mng_ptr->precise_point_num_;
    double precursor_mass = prsm_ptr->getOriPrecMass();
    std::string str = str_util::fixedToString(precursor_mass, pos);
    xml_doc->addElement(ms_header_element, "precursor_mono_mass", str.c_str());

    int precursor_charge = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecCharge();
    str = std::to_string(precursor_charge);
    xml_doc->addElement(ms_header_element, "precursor_charge", str.c_str());

    double precursor_mz = peak_util::compMz(precursor_mass, precursor_charge);
    str = str_util::fixedToString(precursor_mz, pos);
    xml_doc->addElement(ms_header_element, "precursor_mz", str.c_str());

    double feature_inte = prsm_ptr->getFracFeatureInte();
    if (feature_inte > 0) {
      str = str_util::toScientificStr(feature_inte, pos);
      xml_doc->addElement(ms_header_element, "feature_inte", str.c_str());
    }
  }
}

void addMsHeaderBrief(XmlDOMDocument* xml_doc, XmlDOMElement ms_element, 
                 const PrsmPtr &prsm_ptr, const PrsmViewMngPtr &mng_ptr) {
  XmlDOMElement ms_header_element = xml_doc->addElement(ms_element, "ms_header");
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  std::string ms1_ids, ms2_ids;
  std::string ms1_scans, ms2_scans;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    MsHeaderPtr header_ptr = deconv_ms_ptr_vec[i]->getMsHeaderPtr();
    ms1_scans = ms1_scans + std::to_string(header_ptr->getMsOneScan()) + " ";
    ms2_scans = ms2_scans + header_ptr->getScansString() + " ";
  }
  str_util::trim(ms1_scans);
  str_util::trim(ms2_scans);
  xml_doc->addElement(ms_header_element, "ms1_scans", ms1_scans.c_str());
  xml_doc->addElement(ms_header_element, "scans", ms2_scans.c_str());
}

void addMsPeaks(XmlDOMDocument *xml_doc, XmlDOMElement ms_element,
                const PrsmPtr &prsm_ptr, const PrsmViewMngPtr &mng_ptr) {
  // peaks to view
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  ExtendMsPtrVec refine_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  XmlDOMElement peaks = xml_doc->addElement(ms_element, "peaks");
  for (size_t s = 0; s < deconv_ms_ptr_vec.size(); s++) {
    // get ion_pair
    PeakIonPairPtrVec pair_ptrs
        = peak_ion_pair_util::genePeakIonPairs(prsm_ptr->getProteoformPtr(),
                                               refine_ms_ptr_vec[s],
                                               mng_ptr->min_mass_);
    // LOG_DEBUG("pair completed");
    for (size_t i = 0; i < deconv_ms_ptr_vec[s]->size(); i++) {
      XmlDOMElement peak_element = xml_doc->addElement(peaks, "peak");
      std::string str = std::to_string(deconv_ms_ptr_vec[s]->getMsHeaderPtr()->getSpecId());
      xml_doc->addElement(peak_element, "spec_id", str.c_str());

      DeconvPeakPtr peak_ptr = deconv_ms_ptr_vec[s]->getPeakPtr(i);
      str = std::to_string(peak_ptr->getPeakId());
      xml_doc->addElement(peak_element, "peak_id", str.c_str());

      double mass = peak_ptr->getPosition();
      int charge = peak_ptr->getCharge();
      str = str_util::fixedToString(mass, mng_ptr->precise_point_num_);
      xml_doc->addElement(peak_element, "monoisotopic_mass", str.c_str());

      double mz = peak_util::compMz(mass, charge);
      str = str_util::fixedToString(mz, mng_ptr->precise_point_num_);
      xml_doc->addElement(peak_element, "monoisotopic_mz", str.c_str());

      str = str_util::fixedToString(peak_ptr->getIntensity(), mng_ptr->decimal_point_num_);
      xml_doc->addElement(peak_element, "intensity", str.c_str());

      str = std::to_string(charge);
      xml_doc->addElement(peak_element, "charge", str.c_str());

      int spec_id = deconv_ms_ptr_vec[s]->getMsHeaderPtr()->getSpecId();
      PeakIonPairPtrVec selected_pair_ptrs
          = peak_ion_pair_util::getMatchedPairs(pair_ptrs, spec_id, peak_ptr->getPeakId());
      if (selected_pair_ptrs.size() > 0) {
        int match_ions_number = selected_pair_ptrs.size();
        str = std::to_string(match_ions_number);
        xml_doc->addElement(peak_element, "matched_ions_num", str.c_str());

        XmlDOMElement mi_element = xml_doc->addElement(peak_element, "matched_ions");
        for (size_t j = 0; j < selected_pair_ptrs.size(); j++) {
          selected_pair_ptrs[j]->appendTheoPeakToXml(xml_doc, mi_element);
        }
      }
    }
  }
}

XmlDOMElement geneAnnoPrsm(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                  const PrsmPtr &prsm_ptr,
                                  const PrsmViewMngPtr &mng_ptr, bool detail, bool add_ms_peaks) {
  XmlDOMElement prsm_element = xml_doc->addElement(parent, "prsm");
  addPrsmHeader(xml_doc, prsm_element, prsm_ptr, mng_ptr);

  if (detail) {
    XmlDOMElement ms2_element = xml_doc->addElement(prsm_element, "ms");
    addMsHeader(xml_doc, ms2_element, prsm_ptr, mng_ptr);

    if (add_ms_peaks) {
      // add ms peaks
      addMsPeaks(xml_doc, ms2_element, prsm_ptr, mng_ptr);
    }

    // proteoform to view
    anno_proteoform::geneAnnoProteoform(xml_doc, prsm_element, prsm_ptr, mng_ptr);
    LOG_DEBUG("proteoform view completed");
  }
  return prsm_element;
}

XmlDOMElement geneAnnoPrsmBrief(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                  const PrsmPtr &prsm_ptr,
                                  const PrsmViewMngPtr &mng_ptr, bool detail, bool add_ms_peaks) {
  //for prsms.js file in the html folder. Will only contain information needed by TopMSV spectrum id page
  XmlDOMElement prsm_element = xml_doc->addElement(parent, "prsm");
  addPrsmHeaderBrief(xml_doc, prsm_element, prsm_ptr, mng_ptr);

  if (detail) {
    XmlDOMElement ms2_element = xml_doc->addElement(prsm_element, "ms");
    addMsHeaderBrief(xml_doc, ms2_element, prsm_ptr, mng_ptr);


    // proteoform to view
    anno_proteoform::geneAnnoProteoformBrief(xml_doc, prsm_element, prsm_ptr, mng_ptr);
    LOG_DEBUG("proteoform view completed");
  }
  return prsm_element;
}
}
}  // namespace toppic
