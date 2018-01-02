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

#include "base/logger.hpp"
#include "base/ptm_base.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod_base.hpp"
#include "base/prot_mod_util.hpp"
#include "base/trunc_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

namespace prot_mod_util {

bool allowMod(ProtModPtr prot_mod_ptr, const ResiduePtrVec &residues) {
  if (prot_mod_ptr == ProtModBase::getProtModPtr_NONE()) {
    return true;
  } else if (prot_mod_ptr == ProtModBase::getProtModPtr_M_ACETYLATION()) {
    int mod_pos = prot_mod_ptr->getModPos();
    if (mod_pos >= static_cast<int>(residues.size())) {
      // LOG_DEBUG("pos false");
      return false;
    }
    ModPtr mod_ptr = prot_mod_ptr->getModPtr();
    if (residues[mod_pos] != mod_ptr->getOriResiduePtr()) {
      // LOG_DEBUG("mod false");
      return false;
    }
    return true;
  } else {
    // check trunc
    if (!trunc_util::isValidTrunc(prot_mod_ptr->getTruncPtr(), residues)) {
      return false;
    }
    ModPtr mod_ptr = prot_mod_ptr->getModPtr();
    if (mod_ptr != ModBase::getNoneModPtr()) {
      // if NME_acetylation
      int mod_pos = prot_mod_ptr->getModPos();
      if (mod_pos >= static_cast<int>(residues.size())) {
        // LOG_DEBUG("pos false");
        return false;
      }
      if (residues[mod_pos] != mod_ptr->getOriResiduePtr()) {
        // LOG_DEBUG("mod false");
        return false;
      }
    }
    return true;
  }
}

ProtModPtrVec readProtMod(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  ProtModPtrVec mod_ptr_vec;
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = ProtMod::getXmlElementName();
    int mod_num = xml_dom_util::getChildCount(parent, element_name.c_str());
    LOG_DEBUG("mod num " << mod_num);
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element
          = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
      ProtModPtr ptr = ProtModBase::getProtModPtrFromXml(element);
      mod_ptr_vec.push_back(ptr);
    }
  }
  return mod_ptr_vec;
}

ProtModPtr findNME_Acetylation(const ProtModPtrVec &prot_mod_ptrs,
                               const ResiduePtrVec &residues) {
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    PtmPtr ptm_ptr = prot_mod_ptrs[i]->getModPtr()->getModResiduePtr()->getPtmPtr();
    // LOG_DEBUG("ptm ptr " << ptm_ptr->getAbbrName() <<
    //          " equal " << (ptm_ptr == PtmBase::getPtmPtr_Acetylation()));
    if (ptm_ptr == PtmBase::getPtmPtr_Acetylation() &&
        allowMod(prot_mod_ptrs[i], residues)) {
      return prot_mod_ptrs[i];
    }
  }
  return nullptr;
}

} // namespace prot_mod_util

}  // namespace prot
