//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include <string>
#include <vector>

#include "prsm/peak_ion_pair_util.hpp"
#include "prsmview/anno_cleavage.hpp"

namespace toppic {

AnnoCleavage::AnnoCleavage(int pos, const PeakIonPairPtrVec &pairs, 
                           bool exist_n_ion, bool exist_c_ion):
    pos_(pos),
    pairs_(pairs),
    exist_n_ion_(exist_n_ion),
    exist_c_ion_(exist_c_ion) {}

void AnnoCleavage::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("cleavage");
  std::string str = str_util::toString(pos_);
  xml_doc->addElement(element, "position", str.c_str());

  str = str_util::toString(exist_n_ion_);
  xml_doc->addElement(element, "exist_n_ion", str.c_str());

  str = str_util::toString(exist_c_ion_);
  xml_doc->addElement(element, "exist_c_ion", str.c_str());

  xercesc::DOMElement* peaks = xml_doc->createElement("matched_peaks");
  for (size_t i = 0; i < pairs_.size(); i++) {
    pairs_[i]->appendRealPeakToXml(xml_doc, peaks);
  }
  element->appendChild(peaks);
  parent->appendChild(element);
}

AnnoCleavagePtrVec AnnoCleavage::getProteoCleavage(PrsmPtr prsm_ptr, double min_mass) {
  AnnoCleavagePtrVec cleavages;
  ProteoformPtr proteo_ptr = prsm_ptr->getProteoformPtr();
  ExtendMsPtrVec refine_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  int prot_len = proteo_ptr->getFastaSeqPtr()->getAcidPtmPairLen();
  std::vector<bool> n_ion(prot_len + 1, false);
  std::vector<bool> c_ion(prot_len + 1, false);
  PeakIonPairPtrVec2D peak_list(prot_len + 1, PeakIonPairPtrVec(0));

  for (size_t m = 0; m < refine_ms_ptr_vec.size(); m++) {
    PeakIonPairPtrVec pairs
        = peak_ion_pair_util::genePeakIonPairs(proteo_ptr, refine_ms_ptr_vec[m], min_mass);
    for (size_t i = 0; i < pairs.size(); i++) {
      int pos = pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos()+ proteo_ptr->getStartPos();
      // LOG_DEBUG("start pos " << prot_ptr->getStartPos() << " pos " << pos);
      peak_list[pos].push_back(pairs[i]);
      if (pairs[i]->getTheoPeakPtr()->getIonPtr()->getIonTypePtr()->isNTerm()) {
        n_ion[pos] = true;
      } else {
        c_ion[pos] = true;
      }
    }
  }

  for (int i = 0; i < prot_len + 1; i++) {
    AnnoCleavagePtr cleavage
        = std::make_shared<AnnoCleavage>(i, peak_list[i], n_ion[i], c_ion[i]);
    cleavages.push_back(cleavage);
    // LOG_DEBUG("i  " << i << " n ion " << n_ion[i] << " c ion " << c_ion[i]);
  }
  return cleavages;
}


} /* namespace toppic */
