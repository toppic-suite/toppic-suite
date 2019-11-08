//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsmview/anno_proteoform.hpp"
#include "prsmview/anno_prsm.hpp"

namespace toppic {

namespace anno_prsm {

void addPrsmHeader(XmlDOMDocument* xml_doc, xercesc::DOMElement* element,
                   PrsmPtr prsm_ptr, PrsmViewMngPtr mng_ptr) {
  std::string str = str_util::toString(prsm_ptr->getPrsmId());
  xml_doc->addElement(element, "prsm_id", str.c_str());
  if (prsm_ptr->getExtremeValuePtr().get() != nullptr) {
    str = str_util::toString(prsm_ptr->getExtremeValuePtr()->getPValue(), 
                             mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "p_value", str.c_str());
  } else {
    xml_doc->addElement(element, "p_value", "N/A");
  }

  if (prsm_ptr->getExtremeValuePtr().get() != nullptr) {
    str = str_util::toString(prsm_ptr->getExtremeValuePtr()->getEValue(), 
                             mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "e_value", str.c_str());
  } else {
    xml_doc->addElement(element, "e_value", "N/A");
  }

  double fdr = prsm_ptr->getFdr();

  if (fdr >= 0) {
    str = str_util::toString(prsm_ptr->getFdr(), mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "fdr", str.c_str());
  } else {
    xml_doc->addElement(element, "fdr", "N/A");
  }

  str = str_util::toString(static_cast<int>(prsm_ptr->getMatchFragNum()));

  xml_doc->addElement(element, "matched_fragment_number", str.c_str());

  str = str_util::toString(static_cast<int>(prsm_ptr->getMatchPeakNum()));

  xml_doc->addElement(element, "matched_peak_number", str.c_str());
}

void addMsHeader(XmlDOMDocument* xml_doc, xercesc::DOMElement* ms_element, 
                 PrsmPtr prsm_ptr, PrsmViewMngPtr mng_ptr) {
  xercesc::DOMElement* ms_header_element = xml_doc->createElement("ms_header");
  ms_element->appendChild(ms_header_element);
  xml_doc->addElement(ms_header_element, "ms/spectrum_file_name", prsm_ptr->getFileName().c_str());
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  std::string ms1_ids, ms2_ids;
  std::string ms1_scans, ms2_scans;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    ms1_ids = ms1_ids + str_util::toString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getMsOneId()) + " ";
    ms1_scans = ms1_scans + str_util::toString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getMsOneScan()) + " ";
    ms2_ids = ms2_ids + str_util::toString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getId()) + " ";
    ms2_scans = ms2_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
  }
  str_util::trim(ms1_ids);
  str_util::trim(ms1_scans);
  str_util::trim(ms2_ids);
  str_util::trim(ms2_scans);
  xml_doc->addElement(ms_header_element, "ms1_ids", ms1_ids.c_str());
  xml_doc->addElement(ms_header_element, "ms1_scans", ms1_scans.c_str());
  xml_doc->addElement(ms_header_element, "ids", ms2_ids.c_str());
  xml_doc->addElement(ms_header_element, "scans", ms2_scans.c_str());

  int pos = mng_ptr->precise_point_num_;

  double precursor_mass = prsm_ptr->getOriPrecMass();
  std::string str = str_util::toString(precursor_mass, pos);
  xml_doc->addElement(ms_header_element, "precursor_mono_mass", str.c_str());

  int precursor_charge = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge();
  str = str_util::toString(precursor_charge);
  xml_doc->addElement(ms_header_element, "precursor_charge", str.c_str());

  double precursor_mz = Peak::compMz(precursor_mass, precursor_charge);
  str = str_util::toString(precursor_mz, pos);
  xml_doc->addElement(ms_header_element, "precursor_mz", str.c_str());

  double precursor_inte = prsm_ptr->getPrecFeatureInte();
  if (precursor_inte > 0) {
    str = str_util::toScientificStr(precursor_inte, pos);
    xml_doc->addElement(ms_header_element, "precursor_inte", str.c_str());
  }
}

void addMsPeaks(XmlDOMDocument *xml_doc, xercesc::DOMElement* ms_element,
                PrsmPtr prsm_ptr, PrsmViewMngPtr mng_ptr) {
  // peaks to view
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  ExtendMsPtrVec refine_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  xercesc::DOMElement* peaks = xml_doc->createElement("peaks");
  ms_element->appendChild(peaks);
  for (size_t s = 0; s < deconv_ms_ptr_vec.size(); s++) {
    // get ion_pair
    PeakIonPairPtrVec pair_ptrs
        = peak_ion_pair_util::genePeakIonPairs(prsm_ptr->getProteoformPtr(),
                                               refine_ms_ptr_vec[s],
                                               mng_ptr->min_mass_);
    // LOG_DEBUG("pair completed");
    for (size_t i = 0; i < deconv_ms_ptr_vec[s]->size(); i++) {
      xercesc::DOMElement* peak_element = xml_doc->createElement("peak");
      peaks->appendChild(peak_element);
      std::string str = str_util::toString(deconv_ms_ptr_vec[s]->getMsHeaderPtr()->getId());
      xml_doc->addElement(peak_element, "ms/spec_id", str.c_str());

      DeconvPeakPtr peak_ptr = deconv_ms_ptr_vec[s]->getPeakPtr(i);
      str = str_util::toString(peak_ptr->getId());
      xml_doc->addElement(peak_element, "peak_id", str.c_str());

      double mass = peak_ptr->getPosition();
      int charge = peak_ptr->getCharge();
      str = str_util::toString(mass, mng_ptr->precise_point_num_);
      xml_doc->addElement(peak_element, "monoisotopic_mass", str.c_str());

      double mz = Peak::compMz(mass, charge);
      str = str_util::toString(mz, mng_ptr->precise_point_num_);
      xml_doc->addElement(peak_element, "monoisotopic_mz", str.c_str());

      str = str_util::toString(peak_ptr->getIntensity(), mng_ptr->decimal_point_num_);
      xml_doc->addElement(peak_element, "intensity", str.c_str());

      str = str_util::toString(charge);
      xml_doc->addElement(peak_element, "charge", str.c_str());

      int spec_id = deconv_ms_ptr_vec[s]->getMsHeaderPtr()->getId();
      PeakIonPairPtrVec selected_pair_ptrs
          = peak_ion_pair_util::getMatchedPairs(pair_ptrs, spec_id, peak_ptr->getId());
      if (selected_pair_ptrs.size() > 0) {
        int match_ions_number = selected_pair_ptrs.size();
        str = str_util::toString(match_ions_number);
        xml_doc->addElement(peak_element, "matched_ions_num", str.c_str());

        xercesc::DOMElement* mi_element = xml_doc->createElement("matched_ions");
        peak_element->appendChild(mi_element);
        for (size_t j = 0; j < selected_pair_ptrs.size(); j++) {
          selected_pair_ptrs[j]->appendTheoPeakToXml(xml_doc, mi_element);
        }
      }
    }
  }
}

xercesc::DOMElement* geneAnnoPrsm(XmlDOMDocument* xml_doc, PrsmPtr prsm_ptr,
                                  PrsmViewMngPtr mng_ptr, bool detail, bool add_ms_peaks) {
  xercesc::DOMElement* prsm_element = xml_doc->createElement("prsm");
  addPrsmHeader(xml_doc, prsm_element, prsm_ptr, mng_ptr);

  if (detail) {
    xercesc::DOMElement* ms2_element = xml_doc->createElement("ms");
    addMsHeader(xml_doc, ms2_element, prsm_ptr, mng_ptr);

    if (add_ms_peaks) {
      // add ms peaks
      addMsPeaks(xml_doc, ms2_element, prsm_ptr, mng_ptr);
    }

    prsm_element->appendChild(ms2_element);

    // proteoform to view
    xercesc::DOMElement* prot_element 
        = anno_proteoform::geneAnnoProteoform(xml_doc, prsm_ptr, mng_ptr);
    prsm_element->appendChild(prot_element);
    LOG_DEBUG("proteoform view completed");
  }
  return prsm_element;
}

}
}  // namespace toppic
