//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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
#include <cmath>
#include <fstream>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/amino_acid_base.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/mod_base.hpp"
#include "common/base/mod_util.hpp"

namespace toppic {

namespace mod_util {

ModPtrVec readModXml(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  ModPtrVec mod_ptr_vec;
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    XmlDOMElement* parent = doc.getDocumentElement();
    std::string element_name = Mod::getXmlElementName();
    int mod_num = xml_dom_util::getChildCount(parent, element_name.c_str());
    LOG_DEBUG("mod num " << mod_num);
    for (int i = 0; i < mod_num; i++) {
      XmlDOMElement* element
          = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
      ModPtr ptr = ModBase::getModPtrFromXml(element);
      mod_ptr_vec.push_back(ptr);
    }
  }
  return mod_ptr_vec;
}

std::vector<std::vector<std::string>> readModTxtForTsv(const std::string &file_name) {
  //parse mod text file to get mass and mod name to be written in a tsv file
  std::vector<std::vector<std::string>> mod_data;
  std::ifstream infile(file_name.c_str());
  if (!infile.is_open()) {
    LOG_ERROR("PTM file " << file_name <<  "cannot be opened!");
    exit(EXIT_FAILURE);
  }
  std::string line;
  while (std::getline(infile, line)) {
    if (line.size() == 0) continue;
    if (line[0] == '#') continue;
    line = str_util::rmComment(line);
    if (line == "") continue;
    try {
      std::vector<std::string> l = str_util::split(line, ",");
      if (l.size() != 5) throw "The number of commas is not 4.";
      
      std::string mod_name;
      std::string mass;
      
      try {
        mod_name = l[0];
      }
      catch(std::invalid_argument& e){
        throw "Error in PTM name.";
      }
      try {
        mass = l[1];
      }
      catch(std::invalid_argument& e){
        throw "Error in PTM mass.";
      }
      // amino acid
      if (l[2] == "*") l[2] = "ARNDCEQGHILKMFPSTWYV";
      std::vector<std::string> single_mod{mod_name, mass, l[2]};
      mod_data.push_back(single_mod);
    } catch (char const* e) {
      LOG_ERROR("Errors in the PTM file: " << file_name);
      LOG_ERROR("Please check the line:" << line);
      LOG_ERROR("Error message: " << e);
      exit(EXIT_FAILURE);
    }
  }
  return mod_data;
}

ModPtrVec2D readModTxt(const std::string &file_name) {
  LOG_DEBUG("mod txt file " << file_name);
  std::vector<ModPtrVec> mod_ptr_vec2d(3);
  std::ifstream infile(file_name.c_str());
  if (!infile.is_open()) {
    LOG_ERROR("Variable PTM file " << file_name <<  "cannot be opened!");
    exit(EXIT_FAILURE);
  }
  std::string line;
  while (std::getline(infile, line)) {
    if (line.size() == 0) continue;
    if (line[0] == '#') continue;
    line = str_util::rmComment(line);
    if (line == "") continue;
    try {
      std::vector<std::string> l = str_util::split(line, ",");
      if (l.size() != 5) throw "The number of commas is not 4.";

      if (l[2] == "*" && l[3] == "any") throw "* and any cannot be used for the same modification.";
      if (l[2] == "*") l[2] = "ARNDCEQGHILKMFPSTWYV";

      double mass;
      int unimod_id;
      try {
        mass = std::stod(l[1]);
      }
      catch(std::invalid_argument& e){
        throw "Error in PTM mass.";
      }
      try {
        unimod_id = std::stoi(l[4]);
      }
      catch(std::invalid_argument& e){
        throw "Error in UniModId";
      }
      PtmPtr p = std::make_shared<Ptm>(l[0], l[0], mass, unimod_id);
      p = PtmBase::getPtmPtr(p);

      for (size_t i = 0; i < l[2].length(); i++) {
        std::string aa = l[2].substr(i, 1);
        AminoAcidPtr a = AminoAcidBase::getAminoAcidPtrByOneLetter(aa);
        if (a == nullptr) throw "Error in the list of residues.";
        ResiduePtr ori_residue_ptr 
            = ResidueBase::getBaseResiduePtr(
                std::make_shared<Residue>(a, PtmBase::getEmptyPtmPtr()));
        ResiduePtr mod_residue_ptr 
            = ResidueBase::getBaseResiduePtr(std::make_shared<Residue>(a, p));
        ModPtr m 
            = ModBase::getBaseModPtr(std::make_shared<Mod>(ori_residue_ptr, mod_residue_ptr));
        if (l[3] == "N-term") {
          mod_ptr_vec2d[0].push_back(m);
        } else if (l[3] == "C-term") {
          mod_ptr_vec2d[1].push_back(m);
        } else if (l[3] == "any") {
          mod_ptr_vec2d[0].push_back(m);
          mod_ptr_vec2d[1].push_back(m);
          mod_ptr_vec2d[2].push_back(m);
        } else {
          throw "Error in Position.";
        }
      }
    } catch (char const* e) {
      LOG_ERROR("Errors in the Variable PTM file: " << file_name);
      LOG_ERROR("Please check the line: " << line);
      LOG_ERROR("Error message: " << e);
      exit(EXIT_FAILURE);
    }
  }
  infile.close();
  return mod_ptr_vec2d;
}

ModPtrVec readAnywhereModTxt (const std::string &file_name) {
  ModPtrVec2D mod_vec_2d = readModTxt(file_name);
  // the vector at index 2 is the anywhere modifications
  return mod_vec_2d[2];
}

std::vector<double> readModTxtToShiftList(const std::string &file_name) {
  ModPtrVec2D mod_ptr_list_2d = readModTxt(file_name);
  // consider only modifications with the property anywhere
  ModPtrVec mod_ptr_list = mod_ptr_list_2d[2];
  std::vector<double> shift_list; 
  for (size_t i = 0; i < mod_ptr_list.size(); i++) {
    double shift = mod_ptr_list[i]->getShift();
    // if the shift is not in the list
    if (std::find(shift_list.begin(), shift_list.end(), shift) != shift_list.end()) {
      shift_list.push_back(shift);
    }
  }
  return shift_list;
}

ModPtrVec geneFixedModList(const std::string &str) {
  try {
    if (str == "" || str == "C57" || str == "C58") {
      ModPtrVec mod_ptr_vec;
      if (str == "C57") {
        mod_ptr_vec.push_back(ModBase::getC57ModPtr());
      } else if (str == "C58") {
        mod_ptr_vec.push_back(ModBase::getC58ModPtr());
      }
      return mod_ptr_vec;
    } else {
      return readModTxt(str)[2];
    }
  } catch (const char* e) {
    LOG_ERROR("[Exception]" << e);
    exit(EXIT_FAILURE);
  }
}

ResiduePtrVec geneResidueListWithMod(const ResiduePtrVec & residue_list,
                                     const ModPtrVec & fix_mod_list) {
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

std::vector<double> getModMassVec(const ModPtrVec & var_mod_list) {
  std::vector<double> mod_mass_vec;
  for (size_t i = 0; i < var_mod_list.size(); i++) {
    mod_mass_vec.push_back(var_mod_list[i]->getShift());
  }
  std::sort(mod_mass_vec.begin(), mod_mass_vec.end());
  auto last = std::unique(mod_mass_vec.begin(), mod_mass_vec.end(),
                          [] (double a, double b) {
                            return std::abs(a - b) <= std::pow(10, -4);
                          });
  mod_mass_vec.erase(last, mod_mass_vec.end());
  std::sort(mod_mass_vec.begin(), mod_mass_vec.end());
  return mod_mass_vec;
}

} // namespace mod_util

}  // namespace toppic

