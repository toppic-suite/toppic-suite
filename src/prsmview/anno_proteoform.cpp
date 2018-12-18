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


#include <limits>
#include <set>
#include <string>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include "seq/mass_shift_type.hpp"
#include "base/residue_util.hpp"
#include "base/xml_dom_util.hpp"
#include "base/proteoform_factory.hpp"
#include "spec/peak.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsmview/anno_cleavage.hpp"
#include "prsmview/anno_residue.hpp"
#include "prsmview/anno_ptm.hpp"
#include "prsmview/anno_segment.hpp"
#include "prsmview/prsm_view_mng.hpp"

namespace toppic {

void addSummary(XmlDOMDocument* xml_doc, xercesc::DOMElement *prot_element,
                ProteoformPtr proteoform_ptr, PrsmViewMngPtr mng_ptr) {
  std::string str = string_util::convertToString(proteoform_ptr->getProtId());
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());

  str = string_util::convertToString(proteoform_ptr->getProteoClusterId());
  xml_doc->addElement(prot_element, "proteoform_id", str.c_str());

  str = proteoform_ptr->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());

  str = proteoform_ptr->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());

  double mass = proteoform_ptr->getMass();
  str = string_util::convertToString(mass, mng_ptr->precise_point_num_);
  xml_doc->addElement(prot_element, "proteoform_mass", str.c_str());

  str = string_util::convertToString(proteoform_ptr->getProtModPtr()->isAcetylation());
  xml_doc->addElement(prot_element, "n_acetylation", str.c_str());

  int unexpected_change_number = proteoform_ptr->getMassShiftNum(MassShiftType::UNEXPECTED);
  str = string_util::convertToString(unexpected_change_number);
  xml_doc->addElement(prot_element, "unexpected_change_number", str.c_str());
}

void addAnnoHeader(XmlDOMDocument* xml_doc, xercesc::DOMElement *anno_element,
                   ProteoformPtr proteoform_ptr) {
  std::string str = string_util::convertToString(proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairLen());
  xml_doc->addElement(anno_element, "protein_length", str.c_str());

  str = string_util::convertToString(proteoform_ptr->getStartPos());
  xml_doc->addElement(anno_element, "first_residue_position", str.c_str());

  str = string_util::convertToString(proteoform_ptr->getEndPos());
  xml_doc->addElement(anno_element, "last_residue_position", str.c_str());
}

AnnoResiduePtrVec getAnnoResidues(ProteoformPtr proteoform_ptr, PrsmViewMngPtr mng_ptr) {
  AnnoResiduePtrVec res_ptrs;
  StringPairVec acid_ptm_pairs = proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairVec();
  ModPtrVec fix_mod_list = mng_ptr->prsm_para_ptr_->getFixModPtrVec();
  ResiduePtrVec fasta_residues = residue_util::convertStrToResiduePtrVec(acid_ptm_pairs, fix_mod_list);

  int prot_len = fasta_residues.size();
  for (int i = 0; i < prot_len; i++) {
    res_ptrs.push_back(std::make_shared<AnnoResidue>(fasta_residues[i], i));
  }
  return res_ptrs;
}

void addTruncationAnno(ProteoformPtr proteoform_ptr, AnnoCleavagePtrVec &cleavage_ptrs,
                       AnnoResiduePtrVec res_ptrs) {
  // add information for N-terminal truncation
  int start_pos = proteoform_ptr->getStartPos();
  for (int i =0; i < start_pos; i++) {
    cleavage_ptrs[i]->setType(CLEAVAGE_TYPE_N_TRUNCATION);
    res_ptrs[i]->setType(ANNO_RESIDUE_TYPE_N_TRUNCATION);
  }

  if (start_pos > 0) {
    cleavage_ptrs[start_pos]->setType(CLEAVAGE_TYPE_SEQ_START);
  }
  // LOG_DEBUG("n-trunc completed");

  // add information for C-terminal truncation
  int prot_len = proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairLen();
  int end_pos = proteoform_ptr->getEndPos();
  if (end_pos < prot_len - 1) {
    cleavage_ptrs[end_pos + 1]->setType(CLEAVAGE_TYPE_SEQ_END);
  }

  for (int i = end_pos + 1; i < prot_len; i++) {
    cleavage_ptrs[i+1]->setType(CLEAVAGE_TYPE_C_TRUNCATION);
    res_ptrs[i]->setType(ANNO_RESIDUE_TYPE_C_TRUNCATION);
  }
}

// addKnownPtms for fixed and protein variable PTMs
void addKnownPtms(ProteoformPtr proteoform_ptr, MassShiftPtrVec & shift_ptrs,
                  AnnoPtmPtrVec &ptm_ptrs, AnnoResiduePtrVec &res_ptrs) {
  std::sort(shift_ptrs.begin(), shift_ptrs.end(), MassShift::cmpPosInc);
  int start_pos = proteoform_ptr->getStartPos();
  for (size_t i = 0; i < shift_ptrs.size(); i++) {
    PtmPtr ptm_ptr = shift_ptrs[i]->getChangePtrVec()[0]->getModPtr()->getModResiduePtr()->getPtmPtr();
    if (ptm_ptr->getName() == "No PTM") continue;
    int left_db_bp = shift_ptrs[i]->getLeftBpPos() + start_pos;
    res_ptrs[left_db_bp]->setType(ANNO_RESIDUE_TYPE_KNOWN_CHANGE);
    AnnoPtmPtr existing_ptr
        = AnnoPtm::findPtm(ptm_ptrs, ptm_ptr, shift_ptrs[i]->getTypePtr());
    if (existing_ptr == nullptr) {
      existing_ptr = std::make_shared<AnnoPtm>(ptm_ptr, shift_ptrs[i]->getTypePtr());
      ptm_ptrs.push_back(existing_ptr);
    }

    StringPairVec acid_ptm_pairs = proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairVec();
    std::string acid_letter = acid_ptm_pairs[left_db_bp].first;
    existing_ptr->addOccurence(left_db_bp, acid_letter);
  }
}

void addInsertion(int left_db_bp, int right_db_bp, MassShiftPtr shift_ptr, int color,
                  int & last_right, AnnoSegmentPtrVec & segment_ptrs,
                  AnnoCleavagePtrVec & cleavage_ptrs) {
  int this_left = left_db_bp * 2;
  if (this_left > last_right + 1) {
    AnnoSegmentPtr segment_ptr
        = std::make_shared<AnnoSegment>("EMPTY", last_right + 1, this_left - 1, "", -1);
    segment_ptrs.push_back(segment_ptr);
  }

  int this_right = right_db_bp * 2;
  AnnoSegmentPtr segment_ptr = std::make_shared<AnnoSegment>("SHIFT", this_left , this_right,
                                                             shift_ptr->getSeqStr(), color);
  segment_ptrs.push_back(segment_ptr);
  last_right = this_right;
  cleavage_ptrs[left_db_bp]->setUnexpectedChange(true);
  cleavage_ptrs[left_db_bp]->setUnexpectedChangeColor(color);;
}

void addMod(ProteoformPtr proteoform_ptr, int left_db_bp, int right_db_bp,
            MassShiftPtr shift_ptr, int color, int & last_right,
            AnnoSegmentPtrVec & segment_ptrs, AnnoCleavagePtrVec & cleavage_ptrs,
            AnnoResiduePtrVec & res_ptrs) {
  int this_left = left_db_bp * 2 + 1;
  if (this_left > last_right + 1) {
    AnnoSegmentPtr segment_ptr = std::make_shared<AnnoSegment>("EMPTY", last_right + 1,
                                                               this_left - 1, "", -1);
    segment_ptrs.push_back(segment_ptr);
  }
  int this_right = right_db_bp * 2 - 1;
  AnnoSegmentPtr segment_ptr = std::make_shared<AnnoSegment>("SHIFT", this_left , this_right,
                                                             shift_ptr->getSeqStr(), color);
  std::string anno = "";

  StringPairVec acid_ptm_pairs = proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairVec();
  for (int j = left_db_bp; j < right_db_bp; j++) {
    std::string fasta_seq = proteoform_ptr->getFastaSeqPtr()->getRawSeq();
    std::string acid_letter = acid_ptm_pairs[j].first;
    segment_ptr->addOccurence(j, acid_letter);
  }

  LocalAnnoPtr local_anno = shift_ptr->getChangePtr(0)->getLocalAnno();

  if (local_anno != nullptr) {
    for (size_t j = 0; j < local_anno->getScrVec().size(); j++) {
      segment_ptr->addScr(local_anno->getScrVec()[j]);
    }
    if (local_anno->getPtmPtr() != nullptr) {
      segment_ptr->setPtmPtr(local_anno->getPtmPtr());
    }
    anno = segment_ptr->getResidueAnno();
  }

  segment_ptr->setMassShiftType(shift_ptr->getTypePtr());

  segment_ptrs.push_back(segment_ptr);
  last_right = this_right;

  for (int j = left_db_bp; j < right_db_bp - 1; j++) {
    res_ptrs[j]->setUnexpectedChange(true);
    res_ptrs[j]->setUnexpectedChangeColor(color);
    res_ptrs[j]->setAnno(anno);
    cleavage_ptrs[j + 1]->setUnexpectedChange(true);
    cleavage_ptrs[j + 1]->setUnexpectedChangeColor(color);;
  }
  res_ptrs[right_db_bp - 1]->setUnexpectedChange(true);
  res_ptrs[right_db_bp - 1]->setUnexpectedChangeColor(color);;
  res_ptrs[right_db_bp - 1]->setAnno(anno);
}

AnnoSegmentPtrVec getSegments(ProteoformPtr proteoform_ptr, MassShiftPtrVec & shift_ptrs,
                              AnnoCleavagePtrVec &cleavage_ptrs, AnnoResiduePtrVec &res_ptrs) {
  AnnoSegmentPtrVec segment_ptrs;
  int unexpected_shift_color = 0;
  int last_right = -1;
  int start_pos = proteoform_ptr->getStartPos();
  for (size_t i = 0; i < shift_ptrs.size(); i++) {
    int left_db_bp = shift_ptrs[i]->getLeftBpPos() + start_pos;
    int right_db_bp = shift_ptrs[i]->getRightBpPos() + start_pos;
    if (left_db_bp == right_db_bp) {
      // if the shift does not cover any amino acids, start from the cleavage
      addInsertion(left_db_bp, right_db_bp, shift_ptrs[i], unexpected_shift_color,
                   last_right, segment_ptrs, cleavage_ptrs);
    } else {
      addMod(proteoform_ptr, left_db_bp, right_db_bp, shift_ptrs[i], unexpected_shift_color,
             last_right, segment_ptrs, cleavage_ptrs, res_ptrs);
    }

    unexpected_shift_color = (unexpected_shift_color + 1) % 2;
  }

  // last annochange
  AnnoSegmentPtr segment_ptr = std::make_shared<AnnoSegment>("EMPTY", last_right + 1,
                                                             std::numeric_limits<int>::max(),
                                                             "", -1);
  segment_ptrs.push_back(segment_ptr);
  return segment_ptrs;
}

AnnoSegmentPtrVec removeEmptySegment(AnnoSegmentPtrVec &segment_ptrs) {
  AnnoSegmentPtr non_empty_ptr;
  AnnoSegmentPtrVec selected_segment_ptrs;
  for (size_t i = 0; i < segment_ptrs.size(); i++) {
    AnnoSegmentPtr cur_segment_ptr = segment_ptrs[i];
    if (cur_segment_ptr->getType() != "EMPTY") {
      selected_segment_ptrs.push_back(cur_segment_ptr);
      non_empty_ptr = cur_segment_ptr;
    } else {
      if (non_empty_ptr == nullptr) {
        // first empty segment is kept
        selected_segment_ptrs.push_back(cur_segment_ptr);
      } else {
        non_empty_ptr->setRightPos(cur_segment_ptr->getRightPos());
      }
    }
  }
  return selected_segment_ptrs;
}

xercesc::DOMElement* geneAnnoProteoform(XmlDOMDocument* xml_doc,
                                        PrsmPtr prsm_ptr,
                                        PrsmViewMngPtr mng_ptr) {
  ProteoformPtr proteoform_ptr = prsm_ptr->getProteoformPtr();
  xercesc::DOMElement* prot_element = xml_doc->createElement("annotated_protein");
  addSummary(xml_doc, prot_element, proteoform_ptr, mng_ptr);
  xercesc::DOMElement* anno_element = xml_doc->createElement("annotation");
  prot_element->appendChild(anno_element);
  addAnnoHeader(xml_doc, anno_element, proteoform_ptr);
  // LOG_DEBUG("summary completed");

  AnnoCleavagePtrVec cleavage_ptrs = getProteoCleavage(prsm_ptr, mng_ptr->min_mass_);
  // LOG_DEBUG("cleavage completed");

  // obtain residue_ptrs
  AnnoResiduePtrVec res_ptrs = getAnnoResidues(proteoform_ptr, mng_ptr);
  // LOG_DEBUG("residue completed");

  addTruncationAnno(proteoform_ptr, cleavage_ptrs, res_ptrs);

  AnnoPtmPtrVec anno_ptm_ptrs;

  MassShiftPtrVec shift_ptrs = proteoform_ptr->getMassShiftPtrVec(MassShiftType::FIXED);
  addKnownPtms(proteoform_ptr, shift_ptrs, anno_ptm_ptrs, res_ptrs);

  shift_ptrs = proteoform_ptr->getMassShiftPtrVec(MassShiftType::PROTEIN_VARIABLE);
  addKnownPtms(proteoform_ptr, shift_ptrs, anno_ptm_ptrs, res_ptrs);

  MassShiftPtrVec var_shift_ptrs = proteoform_ptr->getMassShiftPtrVec(MassShiftType::VARIABLE);

  shift_ptrs = proteoform_ptr->getMassShiftPtrVec(MassShiftType::UNEXPECTED);
  shift_ptrs.insert(shift_ptrs.end(), var_shift_ptrs.begin(), var_shift_ptrs.end());
  std::sort(shift_ptrs.begin(), shift_ptrs.end(), MassShift::cmpPosInc);

  AnnoSegmentPtrVec segment_ptrs = getSegments(proteoform_ptr, shift_ptrs,
                                               cleavage_ptrs, res_ptrs);

  // remove EMPTY_CHANGES
  // this step may be removed
  AnnoSegmentPtrVec selected_segment_ptrs = removeEmptySegment(segment_ptrs);
  // LOG_DEBUG("unexpected completed");

  for (size_t i = 0; i < res_ptrs.size(); i++) {
    res_ptrs[i]->appendViewXml(xml_doc, anno_element);
  }

  for (size_t i = 0; i < cleavage_ptrs.size(); i++) {
    cleavage_ptrs[i]->appendXml(xml_doc, anno_element);
  }

  for (size_t i = 0; i < selected_segment_ptrs.size(); i++) {
    selected_segment_ptrs[i]->appendXml(xml_doc, anno_element);
  }

  for (size_t i = 0; i < anno_ptm_ptrs.size(); i++) {
    anno_ptm_ptrs[i]->appendXml(xml_doc, anno_element);
  }

  return prot_element;
}

}  // namespace toppic
