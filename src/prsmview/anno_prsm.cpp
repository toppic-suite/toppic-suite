//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include <set>
#include <boost/algorithm/string.hpp>

#include "base/residue_util.hpp"
#include "base/xml_dom_util.hpp"
#include "spec/peak.hpp"
#include "prsm/peak_ion_pair_factory.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsmview/anno_residue.hpp"
#include "prsmview/anno_prsm.hpp"
#include "prsmview/anno_proteoform.hpp"

namespace prot{

void addPrsmHeader(XmlDOMDocument* xml_doc, xercesc::DOMElement* element, 
                   PrsmPtr prsm_ptr, PrsmViewMngPtr mng_ptr) {
  std::string str = string_util::convertToString(prsm_ptr->getPrsmId());
  xml_doc->addElement(element, "prsm_id", str.c_str());
  if(prsm_ptr->getExtremeValuePtr().get()!=nullptr){
    str=string_util::convertToString(prsm_ptr->getExtremeValuePtr()->getPValue(), mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "p_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "p_value", "N/A");
  }
  if(prsm_ptr->getExtremeValuePtr().get()!=nullptr){
    str=string_util::convertToString(prsm_ptr->getExtremeValuePtr()->getEValue(), mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "e_value", str.c_str());
  }
  else{
    xml_doc->addElement(element, "e_value", "N/A");
  }
  double fdr = prsm_ptr->getFdr();
  if (fdr >= 0) {
    str=string_util::convertToString(prsm_ptr->getFdr(), mng_ptr->decimal_point_num_);
    xml_doc->addElement(element, "fdr", str.c_str());
  }
  else {
    xml_doc->addElement(element, "fdr", "N/A");
  }
  str=string_util::convertToString((int)prsm_ptr->getMatchFragNum());
  xml_doc->addElement(element, "matched_fragment_number", str.c_str());
  str=string_util::convertToString((int)prsm_ptr->getMatchPeakNum());
  xml_doc->addElement(element, "matched_peak_number", str.c_str());
}

void addMsHeader(XmlDOMDocument* xml_doc, xercesc::DOMElement* ms_element, 
                 PrsmPtr prsm_ptr, PrsmViewMngPtr mng_ptr) {
  xercesc::DOMElement* ms_header_element = xml_doc->createElement("ms_header");
  ms_element->appendChild(ms_header_element);
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  std::string spec_ids;
  std::string spec_scans;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    spec_ids = spec_ids + std::to_string(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getId()) + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
  }
  boost::algorithm::trim(spec_ids);
  boost::algorithm::trim(spec_scans);
  xml_doc->addElement(ms_header_element, "ids", spec_ids.c_str());
  xml_doc->addElement(ms_header_element, "scans", spec_scans.c_str());
  int pos = 4;
  double precursor_mass = prsm_ptr->getOriPrecMass();
  std::string str=string_util::convertToString(precursor_mass, pos);
  xml_doc->addElement(ms_header_element, "precursor_mono_mass", str.c_str());
  int precursor_charge = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge();
  str = string_util::convertToString(precursor_charge);
  xml_doc->addElement(ms_header_element, "precursor_charge", str.c_str());
  double precursor_mz = Peak::compMonoMz(precursor_mass, precursor_charge); 
  str = string_util::convertToString(precursor_mz, pos);
  xml_doc->addElement(ms_header_element, "precursor_mz", str.c_str());
  double precursor_inte = prsm_ptr->getPrecFeatureInte();
  if (precursor_inte > 0) {
    str = string_util::convertToScientificStr(precursor_inte, pos);
    xml_doc->addElement(ms_header_element, "precursor_inte", str.c_str());
  }
}

void addMsPeaks(XmlDOMDocument *xml_doc, xercesc::DOMElement* ms_element, 
                PrsmPtr prsm_ptr, PrsmViewMngPtr mng_ptr) {
  //peaks to view
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  ExtendMsPtrVec refine_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  xercesc::DOMElement* peaks = xml_doc->createElement("peaks");
  ms_element->appendChild(peaks);
  for (size_t s = 0; s < deconv_ms_ptr_vec.size(); s++) {
    //get ion_pair
    PeakIonPairPtrVec pair_ptrs = PeakIonPairFactory::genePeakIonPairs(prsm_ptr->getProteoformPtr(), 
                                                                       refine_ms_ptr_vec[s],
                                                                       mng_ptr->min_mass_);
    //LOG_DEBUG("pair completed");
    for(size_t i=0;i< deconv_ms_ptr_vec[s]->size();i++){
      xercesc::DOMElement* peak_element = xml_doc->createElement("peak");
      peaks->appendChild(peak_element);
      std::string str = string_util::convertToString(deconv_ms_ptr_vec[s]->getMsHeaderPtr()->getId());
      xml_doc->addElement(peak_element, "spec_id", str.c_str());
      DeconvPeakPtr peak_ptr = deconv_ms_ptr_vec[s]->getPeakPtr(i);
      str=string_util::convertToString(peak_ptr->getId());
      xml_doc->addElement(peak_element, "peak_id", str.c_str());
      double mass = peak_ptr->getPosition();
      int charge = peak_ptr->getCharge();
      str=string_util::convertToString(mass, mng_ptr->precise_point_num_);
      xml_doc->addElement(peak_element, "monoisotopic_mass", str.c_str());
      double mz = Peak::compMonoMz(mass, charge);
      str=string_util::convertToString(mz, mng_ptr->precise_point_num_);
      xml_doc->addElement(peak_element, "monoisotopic_mz", str.c_str());
      str=string_util::convertToString(peak_ptr->getIntensity(), mng_ptr->decimal_point_num_);
      xml_doc->addElement(peak_element, "intensity", str.c_str());
      str=string_util::convertToString(charge);
      xml_doc->addElement(peak_element, "charge", str.c_str());
      int spec_id = deconv_ms_ptr_vec[s]->getMsHeaderPtr()->getId(); 
      PeakIonPairPtrVec selected_pair_ptrs 
          = PeakIonPairUtil::getMatchedPairs(pair_ptrs, spec_id, peak_ptr->getId());
      if(selected_pair_ptrs.size()>0){
        int match_ions_number = selected_pair_ptrs.size();
        str=string_util::convertToString(match_ions_number);
        xml_doc->addElement(peak_element, "matched_ions_num", str.c_str());
        xercesc::DOMElement* mi_element = xml_doc->createElement("matched_ions");
        peak_element->appendChild(mi_element);
        for(size_t j=0;j< selected_pair_ptrs.size();j++){
          selected_pair_ptrs[j]->appendTheoPeakToXml(xml_doc,mi_element);
        }
      }
    }
  }
}

xercesc::DOMElement* geneAnnoPrsm(XmlDOMDocument* xml_doc,PrsmPtr prsm_ptr, 
                                  PrsmViewMngPtr mng_ptr, bool detail){

  xercesc::DOMElement* prsm_element = xml_doc->createElement("prsm");
  addPrsmHeader(xml_doc, prsm_element, prsm_ptr, mng_ptr);

  if (detail) {
    //add ms
    xercesc::DOMElement* ms_element = xml_doc->createElement("ms");
    addMsHeader(xml_doc, ms_element, prsm_ptr, mng_ptr);
    addMsPeaks(xml_doc, ms_element, prsm_ptr, mng_ptr);
    prsm_element->appendChild(ms_element);

    //proteoform to view
    xercesc::DOMElement* prot_element = geneAnnoProteoform(xml_doc, prsm_ptr, mng_ptr);
    prsm_element->appendChild(prot_element);
    LOG_DEBUG("proteoform view completed");
  }
  return prsm_element;
}
}
