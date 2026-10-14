// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "prsm/prsm_deamid_selector.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "common/base/mod.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/ptm.hpp"
#include "common/base/residue.hpp"
#include "common/util/file_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_parser.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "seq/alter.hpp"
#include "seq/alter_type.hpp"
#include "seq/mass_shift.hpp"

namespace toppic {

namespace prsm_deamid_selector {

namespace {

// Deamidation (Unimod 7): monoisotopic mass and the tolerance with which a
// PTM of the variable PTM file is taken for it
const int kDeamidUnimodId = 7;
const double kDeamidMass = 0.984016;
const double kDeamidMassTolerance = 0.001;

// matched fragment numbers that differ by at least this much decide without
// the precursor mass errors
const double kMinFragNumDiff = 5;

// a spectrum's PrSMs of one protein, start and end position
using FormKey = std::tuple<std::string, int, int>;

struct Counts {
  int group_num = 0;
  int deamid_by_frag_num = 0;
  int deamid_by_prec_error_num = 0;
  int zero_by_frag_num = 0;
  int zero_by_prec_error_num = 0;
  int tie_num = 0;
  int removed_num = 0;
};

FormKey getFormKey(const PrsmStrPtr& prsm_ptr) {
  return std::make_tuple(prsm_ptr->getSeqName(),
                         prsm_ptr->getProteoformStartPos(),
                         prsm_ptr->getProteoformEndPos());
}

bool isZeroPtmPrsm(const PrsmStrPtr& prsm_ptr) {
  return prsm_ptr->getUnexpectedPtmNum() == 0 &&
         prsm_ptr->getVariablePtmNum() == 0;
}

// a PrSM without unexpected shifts whose variable PTMs include one of
// deamid_ptm_names
bool isDeamidPrsm(const PrsmStrPtr& prsm_ptr,
                  const std::set<std::string>& deamid_ptm_names) {
  if (prsm_ptr->getUnexpectedPtmNum() != 0 ||
      prsm_ptr->getVariablePtmNum() == 0) {
    return false;
  }
  std::vector<std::string> str_vec = prsm_ptr->getStrVec();
  std::string prsm_str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
  for (size_t i = 0; i < str_vec.size(); i++) {
    prsm_str += str_vec[i];
  }
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser == nullptr) {
    return false;
  }
  XmlDOMDocument doc(parser->parseStr(prsm_str));
  XmlDOMElement root = doc.getDocumentElement();
  XmlDOMElement form_element =
      xml_dom_util::getChildElement(root, "proteoform", 0);
  std::string shift_name = MassShift::getXmlElementName();
  std::string shift_list_name = shift_name + "_list";
  XmlDOMElement list_element =
      xml_dom_util::getChildElement(form_element, shift_list_name.c_str(), 0);
  int shift_num = xml_dom_util::getChildCount(list_element, shift_name.c_str());
  for (int i = 0; i < shift_num; i++) {
    MassShift shift(xml_dom_util::getChildElement(list_element,
                                                  shift_name.c_str(), i));
    for (const AlterPtr& alter_ptr : shift.getAlterPtrVec()) {
      if (alter_ptr->getTypePtr() != AlterType::VARIABLE ||
          alter_ptr->getModPtr() == nullptr) {
        continue;
      }
      const PtmPtr& ptm_ptr =
          alter_ptr->getModPtr()->getModResiduePtr()->getPtmPtr();
      if (ptm_ptr != nullptr &&
          deamid_ptm_names.count(ptm_ptr->getAbbrName()) > 0) {
        return true;
      }
    }
  }
  return false;
}

// the precursor mass (label-free) moved by the whole isotopes nearest to the
// adjusted precursor mass, which is the proteoform mass in the zero-shift and
// the variable PTM search, minus the adjusted precursor mass
double getPrecError(const PrsmStrPtr& prsm_ptr) {
  double label_mass = 0;
  std::string line =
      prsm_util::getXmlLine(prsm_ptr->getStrVec(), "<n_term_label_mass>");
  if (line != "") {
    label_mass = std::stod(prsm_util::getValueStr(line));
  }
  double prec_mass = prsm_ptr->getOriPrecMass() - label_mass;
  double adjusted_mass = prsm_ptr->getAdjustedPrecMass();
  double isotope_mass = mass_constant::getIsotopeMass();
  double offset = std::round((adjusted_mass - prec_mass) / isotope_mass);
  return prec_mass + offset * isotope_mass - adjusted_mass;
}

// the PrSM of the group with the best E-value (the first one on a tie)
PrsmStrPtr getBestPrsm(const PrsmStrPtrVec& prsm_ptrs,
                       const std::vector<size_t>& idxes) {
  PrsmStrPtr best_ptr = prsm_ptrs[idxes[0]];
  for (size_t i = 1; i < idxes.size(); i++) {
    if (prsm_ptrs[idxes[i]]->getEValue() < best_ptr->getEValue()) {
      best_ptr = prsm_ptrs[idxes[i]];
    }
  }
  return best_ptr;
}

bool keepDeamidPrsm(const PrsmStrPtr& zero_ptr, const PrsmStrPtr& deamid_ptr,
                    Counts& counts) {
  double frag_diff =
      deamid_ptr->getMatchTheoFragNum() - zero_ptr->getMatchTheoFragNum();
  if (frag_diff >= kMinFragNumDiff) {
    counts.deamid_by_frag_num++;
    return true;
  }
  if (-frag_diff >= kMinFragNumDiff) {
    counts.zero_by_frag_num++;
    return false;
  }
  double zero_error = std::abs(getPrecError(zero_ptr));
  double deamid_error = std::abs(getPrecError(deamid_ptr));
  if (deamid_error < zero_error) {
    counts.deamid_by_prec_error_num++;
    return true;
  }
  if (zero_error < deamid_error) {
    counts.zero_by_prec_error_num++;
    return false;
  }
  counts.tie_num++;
  return deamid_ptr->getEValue() <= zero_ptr->getEValue();
}

PrsmStrPtrVec selectPrsms(const PrsmStrPtrVec& prsm_ptrs,
                          const std::set<std::string>& deamid_ptm_names,
                          Counts& counts) {
  std::map<FormKey, std::vector<size_t>> zero_idxes;
  std::map<FormKey, std::vector<size_t>> deamid_idxes;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (isZeroPtmPrsm(prsm_ptrs[i])) {
      zero_idxes[getFormKey(prsm_ptrs[i])].push_back(i);
    } else if (isDeamidPrsm(prsm_ptrs[i], deamid_ptm_names)) {
      deamid_idxes[getFormKey(prsm_ptrs[i])].push_back(i);
    }
  }
  std::vector<bool> removed(prsm_ptrs.size(), false);
  for (const auto& zero_entry : zero_idxes) {
    auto deamid_it = deamid_idxes.find(zero_entry.first);
    if (deamid_it == deamid_idxes.end()) {
      continue;
    }
    counts.group_num++;
    PrsmStrPtr zero_ptr = getBestPrsm(prsm_ptrs, zero_entry.second);
    PrsmStrPtr deamid_ptr = getBestPrsm(prsm_ptrs, deamid_it->second);
    const std::vector<size_t>& loser_idxes =
        keepDeamidPrsm(zero_ptr, deamid_ptr, counts) ? zero_entry.second
                                                     : deamid_it->second;
    for (size_t idx : loser_idxes) {
      removed[idx] = true;
      counts.removed_num++;
    }
  }
  PrsmStrPtrVec result_ptrs;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (!removed[i]) {
      result_ptrs.push_back(prsm_ptrs[i]);
    }
  }
  return result_ptrs;
}

}  // namespace

std::set<std::string> getDeamidPtmNames(const std::string& var_ptm_file_name) {
  std::set<std::string> names;
  ModPtrVec2D mod_ptr_vec2d = mod_util::readModTxt(var_ptm_file_name);
  for (const ModPtrVec& mod_ptrs : mod_ptr_vec2d) {
    for (const ModPtr& mod_ptr : mod_ptrs) {
      const PtmPtr& ptm_ptr = mod_ptr->getModResiduePtr()->getPtmPtr();
      if (ptm_ptr->getUnimodId() == kDeamidUnimodId ||
          std::abs(ptm_ptr->getMonoMass() - kDeamidMass) <=
              kDeamidMassTolerance) {
        names.insert(ptm_ptr->getAbbrName());
      }
    }
  }
  return names;
}

void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext,
             const std::set<std::string>& deamid_ptm_names) {
  std::string base_name = file_util::basename(spec_file_name);
  PrsmReader reader(base_name + "." + input_file_ext);
  PrsmXmlWriter writer(base_name + "." + output_file_ext);
  Counts counts;
  PrsmStrPtr prsm_ptr = reader.readOnePrsmStr();
  while (prsm_ptr != nullptr) {
    int spec_id = prsm_ptr->getSpectrumId();
    PrsmStrPtrVec spec_prsm_ptrs;
    while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
      spec_prsm_ptrs.push_back(prsm_ptr);
      prsm_ptr = reader.readOnePrsmStr();
    }
    PrsmStrPtrVec result_ptrs =
        selectPrsms(spec_prsm_ptrs, deamid_ptm_names, counts);
    for (size_t i = 0; i < result_ptrs.size(); i++) {
      writer.write(result_ptrs[i]);
    }
  }
  reader.close();
  writer.close();
  std::cout << "Deamidation PrSM selection - " << counts.group_num
            << " spectrum proteoform positions with a zero-PTM and a "
               "deamidation PrSM; deamidation kept by matched fragments "
            << counts.deamid_by_frag_num << ", by precursor error "
            << counts.deamid_by_prec_error_num
            << "; zero-PTM kept by matched fragments "
            << counts.zero_by_frag_num << ", by precursor error "
            << counts.zero_by_prec_error_num << "; ties " << counts.tie_num
            << "; " << counts.removed_num << " PrSMs removed." << std::endl;
}

}  // namespace prsm_deamid_selector

}  // namespace toppic
