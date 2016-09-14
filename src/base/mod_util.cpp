// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "base/logger.hpp"
#include "base/mod_base.hpp"
#include "base/mod_util.hpp"
#include "base/ptm_base.hpp"
#include "base/acid_base.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ModPtrVec ModUtil::readModXml(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  ModPtrVec mod_ptr_vec;
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Mod::getXmlElementName();
    int mod_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    LOG_DEBUG("mod num " << mod_num);
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element 
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ModPtr ptr = ModBase::getModPtrFromXml(element);
      mod_ptr_vec.push_back(ptr);
    }
  }
  return mod_ptr_vec;
}

std::vector<ModPtrVec> ModUtil::readModTxt(const std::string &file_name) {
  LOG_DEBUG("mod txt file " << file_name);
  std::vector<ModPtrVec> mod_ptr_vec2d(3);
  std::ifstream infile(file_name.c_str());
  if (!infile.is_open()) {
    std::cerr << "Error: variable PTM file "  
        << file_name <<  "can not be opened" << std::endl;
    exit (EXIT_FAILURE);
  }
  std::string line;
  while(std::getline(infile, line)) {
    if (line[0] == '#') continue;
    line = StringUtil::rmComment(line);
    if (line == "") continue;
    try { 
      std::vector<std::string> l = StringUtil::split(line, ',');
      if (l.size() != 5) throw line;

      if (l[2] == "*" && l[3] == "any") throw line;

      if (l[2] == "*") l[2] = "ARNDCEQGHILKMFPSTWYV";

      PtmPtr p = std::make_shared<Ptm>(l[0], l[0], std::stod(l[1]), std::stoi(l[4]), "", "", "");

      p = PtmBase::getPtmPtr(p);

      for (size_t i = 0; i < l[2].length(); i++) {
        AcidPtr a = AcidBase::getAcidPtrByOneLetter(l[2].substr(i, 1));
        ResiduePtr ori_residue_ptr = ResidueBase::getBaseResiduePtr(std::make_shared<Residue>(a, PtmBase::getEmptyPtmPtr()));
        ResiduePtr mod_residue_ptr = ResidueBase::getBaseResiduePtr(std::make_shared<Residue>(a, p));
        ModPtr m = ModBase::getBaseModPtr(std::make_shared<Mod>(ori_residue_ptr, mod_residue_ptr));
        if (l[3] == "N-term") {
          mod_ptr_vec2d[0].push_back(m);
        } else if (l[3] == "C-term") {
          mod_ptr_vec2d[1].push_back(m);
        } else if (l[3] == "any") {
          mod_ptr_vec2d[0].push_back(m);
          mod_ptr_vec2d[1].push_back(m);
          mod_ptr_vec2d[2].push_back(m);
        } else {
          throw line;
        }
      }
    } catch (char const* e) {
      std::cerr << "Errors in the Variable PTM file: " 
          << file_name << std::endl
          << "Please check the line" << std::endl
          << "\t" << e << std::endl;
      exit (EXIT_FAILURE);
    }
  }
  infile.close();
  return mod_ptr_vec2d;
}

ModPtrVec ModUtil::geneFixedModList(const std::string &str) {
  if (str == "" || str == "C57" || str == "C58") {
    ModPtrVec mod_ptr_vec;
    if (str == "C57") {
      mod_ptr_vec.push_back(ModBase::getC57ModPtr());
    }
    else if (str == "C58") {
      mod_ptr_vec.push_back(ModBase::getC58ModPtr());
    }
    return mod_ptr_vec;
  } else {
    return readModTxt(str)[2];
  }
}

ResiduePtrVec ModUtil::geneResidueListWithMod(ResiduePtrVec residue_list,
                                              ModPtrVec fix_mod_list) {
  ResiduePtrVec result;
  for (size_t i = 0; i < residue_list.size(); i++) {
    bool mod = false;
    for (size_t j = 0; j < fix_mod_list.size(); j++) {
      if (fix_mod_list[j]->getOriResiduePtr() == residue_list[i]) {
        mod = true;
        result.push_back(fix_mod_list[j]->getModResiduePtr());
        break;                 
      }
    }
    if (!mod) {
      result.push_back(residue_list[i]);
    }
  }
  return result;
}

} /* end namespace */

