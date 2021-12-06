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

#include <algorithm>

#include "common/base/ptm_base.hpp"
#include "visual/anno_cleavage.hpp"
#include "visual/anno_ptm.hpp"
#include "visual/anno_mass_shift.hpp"
#include "visual/anno_residue.hpp"
#include "visual/anno_proteoform.hpp"

namespace toppic {

namespace anno_proteoform {

void addSummary(XmlDOMDocument* xml_doc, xercesc::DOMElement *prot_element,
                ProteoformPtr proteoform_ptr, PrsmViewMngPtr mng_ptr) {
  std::string str = str_util::toString(proteoform_ptr->getProtId());
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());

  str = str_util::toString(proteoform_ptr->getProteoClusterId());
  xml_doc->addElement(prot_element, "proteoform_id", str.c_str());

  str = proteoform_ptr->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());

  str = proteoform_ptr->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());

  double mass = proteoform_ptr->getMass();
  str = str_util::fixedToString(mass, mng_ptr->precise_point_num_);
  xml_doc->addElement(prot_element, "proteoform_mass", str.c_str());

  str = str_util::toString(proteoform_ptr->getProtModPtr()->isAcetylation());
  xml_doc->addElement(prot_element, "n_acetylation", str.c_str());

  int unexpected_change_number = proteoform_ptr->getAlterNum(AlterType::UNEXPECTED);
  str = str_util::toString(unexpected_change_number);
  xml_doc->addElement(prot_element, "unexpected_shift_number", str.c_str());
}

void addSummaryBrief(XmlDOMDocument* xml_doc, xercesc::DOMElement *prot_element,
                ProteoformPtr proteoform_ptr, PrsmViewMngPtr mng_ptr) {
  std::string str = proteoform_ptr->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());

  str = proteoform_ptr->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());
}

void addAnnoHeader(XmlDOMDocument* xml_doc, xercesc::DOMElement *anno_element,
                   ProteoformPtr proteoform_ptr) {
  std::string str = str_util::toString(proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairLen());
  xml_doc->addElement(anno_element, "protein_length", str.c_str());

  str = str_util::toString(proteoform_ptr->getStartPos());
  xml_doc->addElement(anno_element, "first_residue_position", str.c_str());

  str = str_util::toString(proteoform_ptr->getEndPos());
  xml_doc->addElement(anno_element, "last_residue_position", str.c_str());

  str = proteoform_ptr->getProteoformMatchSeq();
  xml_doc->addElement(anno_element, "annotated_seq", str.c_str());
}

void addAnnoHeaderBrief(XmlDOMDocument* xml_doc, xercesc::DOMElement *anno_element,
                   ProteoformPtr proteoform_ptr) {
  std::string str = proteoform_ptr->getProteoformMatchSeq();
  xml_doc->addElement(anno_element, "annotated_seq", str.c_str());
}

// addKnownPtms for fixed and protein variable PTMs
void addAnnoPtms(ProteoformPtr proteoform_ptr, MassShiftPtrVec & shift_ptrs,
                 AnnoPtmPtrVec &ptm_ptrs) {
  std::sort(shift_ptrs.begin(), shift_ptrs.end(), MassShift::cmpPosInc);
  int start_pos = proteoform_ptr->getStartPos();
  for (size_t i = 0; i < shift_ptrs.size(); i++) {
    PtmPtr ptm_ptr = shift_ptrs[i]->getAlterPtrVec()[0]->getModPtr()->getModResiduePtr()->getPtmPtr();
    if (PtmBase::isEmptyPtmPtr(ptm_ptr)) {
      continue;
    }
    int left_db_bp = shift_ptrs[i]->getLeftBpPos() + start_pos;
    AnnoPtmPtr existing_ptr
        = AnnoPtm::findPtm(ptm_ptrs, ptm_ptr, shift_ptrs[i]->getTypePtr());
    if (existing_ptr == nullptr) {
      existing_ptr = std::make_shared<AnnoPtm>(ptm_ptr, shift_ptrs[i]->getTypePtr());
      ptm_ptrs.push_back(existing_ptr);
    }

    StringPairVec acid_ptm_pairs = proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairVec();
    std::string acid_letter = acid_ptm_pairs[left_db_bp].first;
    existing_ptr->addOccurence(left_db_bp, left_db_bp + 1, acid_letter);
  }
}

void addAnnoVarPtms(ProteoformPtr proteoform_ptr, MassShiftPtrVec & shift_ptrs,
                    AnnoPtmPtrVec &var_ptm_ptrs) {
  std::sort(shift_ptrs.begin(), shift_ptrs.end(), MassShift::cmpPosInc);
  int start_pos = proteoform_ptr->getStartPos();
  for (size_t i = 0; i < shift_ptrs.size(); i++) {
    AlterPtrVec alter_ptrs = shift_ptrs[i]->getAlterPtrVec(); 
    for (size_t j = 0; j < alter_ptrs.size(); j++) {
      PtmPtr ptm_ptr = alter_ptrs[j]->getModPtr()->getModResiduePtr()->getPtmPtr();
      if (PtmBase::isEmptyPtmPtr(ptm_ptr)) {
        continue;
      }
      int left_db_bp = alter_ptrs[j]->getLeftBpPos() + start_pos;
      int right_db_bp = alter_ptrs[j]->getRightBpPos() + start_pos;
      AnnoPtmPtr existing_ptr
          = AnnoPtm::findPtm(var_ptm_ptrs, ptm_ptr, alter_ptrs[j]->getTypePtr());
      if (existing_ptr == nullptr) {
        existing_ptr = std::make_shared<AnnoPtm>(ptm_ptr, alter_ptrs[j]->getTypePtr());
        var_ptm_ptrs.push_back(existing_ptr);
      }

      StringPairVec acid_ptm_pairs = proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairVec();
      existing_ptr->addOccurence(left_db_bp, right_db_bp, "");
    }
  }
}


void addAnnoMassShifts(ProteoformPtr proteoform_ptr, MassShiftPtrVec & shift_ptrs,
                       AnnoMassShiftPtrVec &anno_shift_ptrs) {
  std::sort(shift_ptrs.begin(), shift_ptrs.end(), MassShift::cmpPosInc);
  int start_pos = proteoform_ptr->getStartPos();
  for (size_t i = 0; i < shift_ptrs.size(); i++) {
    int left_db_bp = shift_ptrs[i]->getLeftBpPos() + start_pos;
    int right_db_bp = shift_ptrs[i]->getRightBpPos() + start_pos;
    double mass_shift = shift_ptrs[i]->getMassShift();
    std::string anno_str = shift_ptrs[i]->getAnnoStr();
    AlterTypePtr type_ptr = shift_ptrs[i]->getTypePtr();
    AnnoMassShiftPtr anno_shift_ptr 
        = std::make_shared<AnnoMassShift>(i, left_db_bp, right_db_bp, mass_shift, 
                                          anno_str, type_ptr);
    anno_shift_ptrs.push_back(anno_shift_ptr);
  }
}

/*
void addInsertion(int left_db_bp, int right_db_bp, MassShiftPtr shift_ptr,
                  AnnoMassShiftPerVec & anno_shift_ptrs) {
  AnnoMassShiftPtr anno_shift_ptr 
      = std::make_shared<AnnoSegment>("SHIFT", left_db_bp, right_db_bp,
                                      shift_ptr->getSeqStr());
  segment_ptrs.push_back(segment_ptr);
  last_right = this_right;
}

AnnoMassShiftPtrVec getMassShifts(ProteoformPtr proteoform_ptr, MassShiftPtrVec & shift_ptrs) {
  AnnoMassShiftVec mass_shift_ptrs;
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
  }
  return mass_shift_ptrs;
}



void addMod(ProteoformPtr proteoform_ptr, int left_db_bp, int right_db_bp,
            MassShiftPtr shift_ptr, AnnoSegmentPtrVec & segment_ptrs) {
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

  segment_ptr->setAlterType(shift_ptr->getTypePtr());

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
*/

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

  // obtain residue_ptrs
  AnnoResiduePtrVec res_ptrs = AnnoResidue::getAnnoResidues(proteoform_ptr);
  for (size_t i = 0; i < res_ptrs.size(); i++) {
    res_ptrs[i]->appendViewXml(xml_doc, anno_element);
  }
  // LOG_DEBUG("residue completed");

  AnnoCleavagePtrVec cleavage_ptrs = AnnoCleavage::getProteoCleavage(prsm_ptr, mng_ptr->min_mass_);
  for (size_t i = 0; i < cleavage_ptrs.size(); i++) {
    cleavage_ptrs[i]->appendXml(xml_doc, anno_element);
  }
  // LOG_DEBUG("cleavage completed");

  // Fixed PTMs
  AnnoPtmPtrVec anno_fixed_ptm_ptrs;
  MassShiftPtrVec fixed_shift_ptrs 
      = proteoform_ptr->getMassShiftPtrVec(AlterType::FIXED);
  addAnnoPtms(proteoform_ptr, fixed_shift_ptrs, anno_fixed_ptm_ptrs);
  for (size_t i = 0; i < anno_fixed_ptm_ptrs.size(); i++) {
    anno_fixed_ptm_ptrs[i]->appendXml(xml_doc, anno_element);
  }

  // protein N-terminal variable PTMs
  AnnoPtmPtrVec anno_prot_var_ptm_ptrs;
  MassShiftPtrVec prot_var_shift_ptrs 
      = proteoform_ptr->getMassShiftPtrVec(AlterType::PROTEIN_VARIABLE);
  addAnnoPtms(proteoform_ptr, prot_var_shift_ptrs, anno_prot_var_ptm_ptrs);
  for (size_t i = 0; i < anno_prot_var_ptm_ptrs.size(); i++) {
    anno_prot_var_ptm_ptrs[i]->appendXml(xml_doc, anno_element);
  }
  // LOG_DEBUG("Fixed PTM and protein vairable PTM finished!)

  // variable PTMs
  AnnoPtmPtrVec anno_var_ptm_ptrs;
  MassShiftPtrVec var_shift_ptrs 
      = proteoform_ptr->getMassShiftPtrVec(AlterType::VARIABLE);
  addAnnoVarPtms(proteoform_ptr, var_shift_ptrs, anno_var_ptm_ptrs); 
  for (size_t i = 0; i < anno_var_ptm_ptrs.size(); i++) {
    anno_var_ptm_ptrs[i]->appendXml(xml_doc, anno_element);
  }

  // Bariable and unexpected mass shifts
  MassShiftPtrVec unexpected_shift_ptrs 
      = proteoform_ptr->getMassShiftPtrVec(AlterType::UNEXPECTED);
  unexpected_shift_ptrs.insert(unexpected_shift_ptrs.end(), 
                               var_shift_ptrs.begin(), var_shift_ptrs.end());
  std::sort(unexpected_shift_ptrs.begin(), unexpected_shift_ptrs.end(), MassShift::cmpPosInc);

  AnnoMassShiftPtrVec anno_mass_shift_ptrs;
  addAnnoMassShifts(proteoform_ptr, unexpected_shift_ptrs, anno_mass_shift_ptrs);

  for (size_t i = 0; i < anno_mass_shift_ptrs.size(); i++) {
    anno_mass_shift_ptrs[i]->appendXml(xml_doc, anno_element);
  }

  /*
  AnnoSegmentPtrVec segment_ptrs = getSegments(proteoform_ptr, shift_ptrs,
                                               cleavage_ptrs, res_ptrs);

  // remove EMPTY_CHANGES
  // this step may be removed
  AnnoSegmentPtrVec selected_segment_ptrs = removeEmptySegment(segment_ptrs);
  // LOG_DEBUG("unexpected completed");


  for (size_t i = 0; i < selected_segment_ptrs.size(); i++) {
    selected_segment_ptrs[i]->appendXml(xml_doc, anno_element);
  }
  */
  return prot_element;
}

xercesc::DOMElement* geneAnnoProteoformBrief(XmlDOMDocument* xml_doc,
                                        PrsmPtr prsm_ptr,
                                        PrsmViewMngPtr mng_ptr) {
  ProteoformPtr proteoform_ptr = prsm_ptr->getProteoformPtr();
  xercesc::DOMElement* prot_element = xml_doc->createElement("annotated_protein");
  addSummaryBrief(xml_doc, prot_element, proteoform_ptr, mng_ptr);
  xercesc::DOMElement* anno_element = xml_doc->createElement("annotation");
  prot_element->appendChild(anno_element);
  addAnnoHeaderBrief(xml_doc, anno_element, proteoform_ptr);
  // LOG_DEBUG("summary completed");

  return prot_element;
}
}

}  // namespace toppic
