#include <set>
#include <boost/algorithm/string.hpp>

#include "base/residue_util.hpp"
#include "base/xml_dom_util.hpp"
#include "base/proteoform_factory.hpp"
#include "spec/peak.hpp"
#include "prsm/peak_ion_pair_factory.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsmview/anno_cleavage.hpp"
#include "prsmview/anno_residue.hpp"
#include "prsmview/prsm_view_mng.hpp"

namespace prot{

void addSummary(XmlDOMDocument* xml_doc, xercesc::DOMElement *prot_element, 
                ProteoformPtr proteoform_ptr, PrsmViewMngPtr mng_ptr) {
  //std::string str=StringUtil::convertToString(proteoform_ptr->getSeqId());
  //xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  std::string str=StringUtil::convertToString(proteoform_ptr->getSpeciesId());
  xml_doc->addElement(prot_element, "proteoform_id", str.c_str());
  str=proteoform_ptr->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  str=proteoform_ptr->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());
  double mass = proteoform_ptr->getMass();
  str=StringUtil::convertToString(mass, mng_ptr->decimal_point_num_);
  xml_doc->addElement(prot_element, "proteoform_mass", str.c_str());
  str=StringUtil::convertToString(proteoform_ptr->getProtModPtr()->isAcetylation());
  xml_doc->addElement(prot_element, "n_acetylation", str.c_str());
  int unexpected_change_number = proteoform_ptr->getChangeNum(ChangeType::UNEXPECTED);
  str=StringUtil::convertToString(unexpected_change_number);
  xml_doc->addElement(prot_element, "unexpected_change_number", str.c_str());

  ChangePtrVec change_ptrs = proteoform_ptr->getChangePtrVec(); 
  std::sort(change_ptrs.begin(),change_ptrs.end(),Change::cmpTypeIncPosInc);
}

void addAnnoHeader(XmlDOMDocument* xml_doc, xercesc::DOMElement *anno_element, 
                ProteoformPtr proteoform_ptr, PrsmViewMngPtr mng_ptr) {

  std::string str=StringUtil::convertToString(proteoform_ptr->getFastaSeqPtr()->getLen());
  xml_doc->addElement(anno_element, "protein_length", str.c_str());

  str=StringUtil::convertToString(proteoform_ptr->getStartPos());
  xml_doc->addElement(anno_element, "first_residue_position", str.c_str());
  str=StringUtil::convertToString(proteoform_ptr->getEndPos());
  xml_doc->addElement(anno_element, "last_residue_position", str.c_str());
  str=StringUtil::convertToString(proteoform_ptr->getProtModPtr()->isAcetylation());
}

AnnoResiduePtrVec getAnnoResidues(ProteoformPtr proteoform_ptr, PrsmViewMngPtr mng_ptr) {
  AnnoResiduePtrVec res_ptrs;
  std::string fasta_seq = proteoform_ptr->getFastaSeqPtr()->getSeq();
  ModPtrVec fix_mod_list = mng_ptr->prsm_para_ptr_->getFixModPtrVec();
  ResiduePtrVec fasta_residues = ResidueUtil::convertStrToResiduePtrVec(fasta_seq,fix_mod_list); 
  int prot_len = proteoform_ptr->getFastaSeqPtr()->getLen();
  for(int i=0;i< prot_len;i++){
    res_ptrs.push_back(AnnoResiduePtr(new AnnoResidue(fasta_residues[i], i)));
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
  //LOG_DEBUG("n-trunc completed");

  // add information for C-terminal truncation
  int prot_len = proteoform_ptr->getFastaSeqPtr()->getLen();
  int end_pos = proteoform_ptr->getEndPos();
  if (end_pos < prot_len - 1) {
    cleavage_ptrs[end_pos + 1]->setType(CLEAVAGE_TYPE_SEQ_END);
  }

  for (int i = end_pos + 1; i < prot_len; i++) {
    cleavage_ptrs[i+1]->setType(CLEAVAGE_TYPE_C_TRUNCATION);
    res_ptrs[i]->setType(ANNO_RESIDUE_TYPE_C_TRUNCATION);
  }
}


xercesc::DOMElement* geneProteinView(XmlDOMDocument* xml_doc,
                                     PrsmPtr prsm_ptr,
                                     PrsmViewMngPtr mng_ptr, double err) {
  ProteoformPtr proteoform_ptr = prsm_ptr->getProteoformPtr();
  xercesc::DOMElement* prot_element = xml_doc->createElement("annotated_protein");
  addSummary(xml_doc, prot_element, proteoform_ptr, mng_ptr);
  xercesc::DOMElement* anno_element = xml_doc->createElement("annotation");
  prot_element->appendChild(anno_element);
  addAnnoHeader(xml_doc, anno_element, proteoform_ptr, mng_ptr);
  //LOG_DEBUG("summary completed");

  AnnoCleavagePtrVec cleavage_ptrs = getProteoCleavage(prsm_ptr, mng_ptr->min_mass_);
  //LOG_DEBUG("cleavage completed");

  // obtain residue_ptrs 
  AnnoResiduePtrVec res_ptrs = getAnnoResidues(proteoform_ptr, mng_ptr);
  //LOG_DEBUG("residue completed");
  
  addTruncationAnno(proteoform_ptr, cleavage_ptrs, res_ptrs);

  /*
  AnnoUnexpectedChangePtrVec unexpected_change_ptrs;
  AnnoExpectedChangePtrVec expected_change_ptrs;
  int unexpected_shift_color = 0;
  int last_right = -1;
  for (size_t i = 0; i < change_ptrs.size(); i++) {
    // if the mass is less than 1 Da
    if (std::abs(change_ptrs[i]->getMassShift()) <= 1 + err)
        continue;
    // add information for known changes 
    int left_db_bp = change_ptrs[i]->getLeftBpPos() + start_pos;
    int right_db_bp = change_ptrs[i]->getRightBpPos() + start_pos;
    double shift = change_ptrs[i]->getMassShift();
    if (change_ptrs[i]->getChangeTypePtr() != ChangeType::UNEXPECTED) { 
      res_ptrs[left_db_bp]->setType(ANNO_RESIDUE_TYPE_KNOWN_CHANGE);
      AnnoExpectedChangePtr existing_ptr 
          = findExpectedChange(expected_change_ptrs, change_ptrs[i]->getChangeTypePtr(), change_ptrs[i]->getModPtr());
      if (existing_ptr == nullptr) {
        existing_ptr = AnnoExpectedChangePtr(new AnnoExpectedChange(change_ptrs[i]->getChangeTypePtr(), 
                                                                    change_ptrs[i]->getModPtr()));
        expected_change_ptrs.push_back(existing_ptr);
      }
      std::string fasta_seq = proteoform_ptr->getFastaSeqPtr()->getSeq();
      std::string acid_letter = fasta_seq.substr(left_db_bp, 1);
      existing_ptr->addOccurence(left_db_bp, acid_letter);
    }
    else {
      if (left_db_bp == right_db_bp) {
        int this_left = left_db_bp * 2;
        if (this_left > last_right + 1) {
          AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(last_right + 1, this_left - 1, 0, -1, "EMPTY"));
          unexpected_change_ptrs.push_back(anno_change_ptr);
        }
        int this_right = right_db_bp * 2;
        AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(this_left , this_right, shift, unexpected_shift_color, "SHIFT"));
        anno_change_ptr->setModPtr(change_ptrs[i]->getModPtr());
        std::string anno_info = "PTM: ";
        if (change_ptrs[i]->getModPtr() == nullptr) {
            anno_info += "Unknown";
        } else {
            anno_info += change_ptrs[i]->getModPtr()->getModResiduePtr()->getPtmPtr()->getAbbrName();
        }
        std::string fasta_seq = proteoform_ptr->getFastaSeqPtr()->getSeq();
        for (int k = left_db_bp; k <= right_db_bp; k++) {
            std::string acid_letter = fasta_seq.substr(k,1);
            anno_change_ptr->addOccurence(k, acid_letter);
            res_ptrs[k]->setPossiblePosColor(1);
            res_ptrs[k]->setAnno(anno_info);
        }

        unexpected_change_ptrs.push_back(anno_change_ptr);
        last_right = this_right;
        cleavage_ptrs[left_db_bp]->setUnexpectedChange(true);
        cleavage_ptrs[left_db_bp]->setUnexpectedChangeColor(unexpected_shift_color);;
      }
      else {
        int this_left = left_db_bp * 2 + 1;
        if (this_left > last_right + 1) {
          AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(last_right + 1, this_left - 1, 0, -1, "EMPTY"));
          unexpected_change_ptrs.push_back(anno_change_ptr);
        }
        int this_right = right_db_bp * 2 - 1;
        AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(this_left, this_right, shift, unexpected_shift_color, "SHIFT"));

        anno_change_ptr->setModPtr(change_ptrs[i]->getModPtr());
        std::string anno_info = "PTM: ";
        if (change_ptrs[i]->getModPtr() != nullptr) {
            anno_info += change_ptrs[i]->getModPtr()->getModResiduePtr()->getPtmPtr()->getName() + "\n";
            std::vector<double> scr = change_ptrs[i]->getScr();
            for (int k = left_db_bp; k < right_db_bp; k++) {
                if (scr[k - left_db_bp] > 0) {
                    std::string acid_letter = proteoform_ptr->getDbResSeqPtr()
                        ->getResiduePtr(k)->getAcidPtr()->getOneLetter();
                    anno_info += "Site: " + acid_letter + std::to_string(k) + " ";
                    anno_info += "Confidence: "
                        + StringUtil::convertToString(scr[k - left_db_bp] * 100, 2) + "%\n";
                }
            }
            for (int k = left_db_bp; k < right_db_bp; k++) {
                if (scr[k - left_db_bp] > 0) {
                    std::string acid_letter = proteoform_ptr->getDbResSeqPtr()
                        ->getResiduePtr(k)->getAcidPtr()->getOneLetter();
                    anno_change_ptr->addOccurence(k, acid_letter);
                    res_ptrs[k]->setPossiblePosColor(1);
                    res_ptrs[k]->setAnno(anno_info);
                }
            }
        } else {
            anno_info += "Unknown\n";
            std::string acid_letter = proteoform_ptr->getDbResSeqPtr()
                ->getResiduePtr(left_db_bp)->getAcidPtr()->getOneLetter();
            anno_change_ptr->addOccurence(left_db_bp, acid_letter);
            anno_info += "Region: " + acid_letter + std::to_string(left_db_bp) + " - ";
            acid_letter = proteoform_ptr->getDbResSeqPtr()->getResiduePtr(
                    right_db_bp - 1)->getAcidPtr()->getOneLetter();
            anno_info += acid_letter + std::to_string(right_db_bp - 1);
            anno_change_ptr->addOccurence(right_db_bp - 1, acid_letter);
            double scr_sum = 0.0;
            std::vector<double> scr = change_ptrs[i]->getScr();
            for (int k = left_db_bp; k < right_db_bp; k++) {
                scr_sum += scr[k - left_db_bp];
            }
            anno_info += " Confindence: " + StringUtil::convertToString(scr_sum * 100, 2) + "%\n";
            for (int k = left_db_bp; k < right_db_bp; k++) {
                res_ptrs[k]->setPossiblePosColor(1);
                res_ptrs[k]->setAnno(anno_info);
            }

        }

        unexpected_change_ptrs.push_back(anno_change_ptr);
        last_right = this_right;

        for (int j = left_db_bp; j < right_db_bp - 1; j++) {
          res_ptrs[j]->setUnexpectedChange(true);
          res_ptrs[j]->setUnexpectedChangeColor(unexpected_shift_color);;
          cleavage_ptrs[j+1]->setUnexpectedChange(true);
          cleavage_ptrs[j+1]->setUnexpectedChangeColor(unexpected_shift_color);;
        }
        res_ptrs[right_db_bp-1]->setUnexpectedChange(true);
        res_ptrs[right_db_bp-1]->setUnexpectedChangeColor(unexpected_shift_color);;
      }
      unexpected_shift_color = (unexpected_shift_color) + 1 % 2;
    }
  }
  // last annochange
  AnnoUnexpectedChangePtr anno_change_ptr(new AnnoUnexpectedChange(last_right + 1, std::numeric_limits<int>::max(), 0, -1, "EMPTY"));
  unexpected_change_ptrs.push_back(anno_change_ptr);
  // remove EMPTY_CHANGES 
  AnnoUnexpectedChangePtr non_empty_ptr;
  AnnoUnexpectedChangePtrVec short_unexpected_change_ptrs;
  for (size_t i = 0; i < unexpected_change_ptrs.size(); i++) {
    AnnoUnexpectedChangePtr cur_change_ptr = unexpected_change_ptrs[i];
    if (cur_change_ptr->getChangeType() != "EMPTY") {
      short_unexpected_change_ptrs.push_back(cur_change_ptr);
      non_empty_ptr = cur_change_ptr;
    }
    else {
      if (non_empty_ptr == nullptr) {
        // first empty segment is kept 
        short_unexpected_change_ptrs.push_back(cur_change_ptr);
      }
      else {
        non_empty_ptr->setRightPos(cur_change_ptr->getRightPos());
      }
    }
  }

  //LOG_DEBUG("unexpected completed");

  for (size_t i = 0; i < res_ptrs.size(); i++) {
    res_ptrs[i]->appendViewXml(xml_doc, anno_element);
  }
  for (size_t i = 0; i < cleavage_ptrs.size(); i++) {
    cleavage_ptrs[i]->appendXml(xml_doc, anno_element);
  }
  for (size_t i = 0; i < short_unexpected_change_ptrs.size(); i++) {
    short_unexpected_change_ptrs[i]->appendXml(xml_doc, anno_element, mng_ptr->decimal_point_num_);
  }
  for (size_t i = 0; i < expected_change_ptrs.size(); i++) {
    expected_change_ptrs[i]->appendXml(xml_doc, anno_element);
  }
  */
  return prot_element;
}

}
