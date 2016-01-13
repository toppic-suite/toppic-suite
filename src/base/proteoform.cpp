#include <sstream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/change_type.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod_base.hpp"
#include "base/string_util.hpp"
#include "base/fasta_index_reader.hpp"
#include "base/proteoform.hpp"
#include "base/proteoform_factory.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Proteoform::Proteoform(FastaSeqPtr fasta_seq_ptr, ProtModPtr prot_mod_ptr, 
                       int start_pos, int end_pos, ResSeqPtr res_seq_ptr, 
                       const ChangePtrVec &change_ptr_vec):
    fasta_seq_ptr_(fasta_seq_ptr),
    prot_mod_ptr_(prot_mod_ptr),
    start_pos_(start_pos),
    end_pos_(end_pos),
    residue_seq_ptr_(res_seq_ptr),
    change_list_(change_ptr_vec) {
      bp_spec_ptr_ = BpSpecPtr(new BpSpec(res_seq_ptr));
      std::sort(change_list_.begin(), change_list_.end(), Change::cmpPosInc);
      species_id_=0;
    }

Proteoform::Proteoform(xercesc::DOMElement* element, FastaIndexReaderPtr reader_ptr,
                       const ModPtrVec &fix_mod_list) {

  std::string seq_element_name = FastaSeq::getXmlElementName();
  xercesc::DOMElement* seq_element= XmlDomUtil::getChildElement(element,seq_element_name.c_str(),0);
  std::string seq_name = FastaSeq::getNameFromXml(seq_element);
  std::string seq_desc = FastaSeq::getDescFromXml(seq_element);

  ProteoformPtr form_ptr = ProteoformFactory::readFastaToProteoformPtr(reader_ptr, seq_name,
                                                                       seq_desc, fix_mod_list);
  parseXml(element, form_ptr);
}

void Proteoform::parseXml(xercesc::DOMElement* element, ProteoformPtr form_ptr) {
  //LOG_DEBUG("start parse proteoform");
  start_pos_ = XmlDomUtil::getIntChildValue(element, "start_pos", 0);
  end_pos_ = XmlDomUtil::getIntChildValue(element, "end_pos", 0);
  species_id_ = XmlDomUtil::getIntChildValue(element, "species_id", 0);
  prot_id_ = XmlDomUtil::getIntChildValue(element, "prot_id", 0);

  //LOG_DEBUG("start parse prot mod");
  std::string pm_element_name = ProtMod::getXmlElementName();
  xercesc::DOMElement* pm_element= XmlDomUtil::getChildElement(element,pm_element_name.c_str(),0);
  prot_mod_ptr_ = ProtModBase::getProtModPtrFromXml(pm_element);

  fasta_seq_ptr_ = form_ptr->getFastaSeqPtr();
  residue_seq_ptr_ = form_ptr->getResSeqPtr()->getSubResidueSeq(start_pos_,end_pos_);

  ModPtr mod_ptr = prot_mod_ptr_->getModPtr();
  if (!ModBase::isNoneModPtr(mod_ptr)) {
    if (residue_seq_ptr_->getLen() >= 1 
        && mod_ptr->getOriResiduePtr() == residue_seq_ptr_->getResiduePtr(0)) {
      ResiduePtrVec residues = residue_seq_ptr_->getResidues();
      residues[0]=mod_ptr->getModResiduePtr();
      residue_seq_ptr_ = ResSeqPtr(new ResidueSeq(residues));

    }
  }

  bp_spec_ptr_= BpSpecPtr(new BpSpec(residue_seq_ptr_));

  //LOG_DEBUG("start parse changes");
  std::string change_element_name = Change::getXmlElementName();

  xercesc::DOMElement* change_list_element= XmlDomUtil::getChildElement(element,"change_list",0);
  int change_len = XmlDomUtil::getChildCount(change_list_element, change_element_name.c_str());

  for(int i=0; i<change_len; i++) {
    xercesc::DOMElement* change_element 
        = XmlDomUtil::getChildElement(change_list_element, change_element_name.c_str(), i);
    change_list_.push_back(ChangePtr(new Change(change_element)));
  }
  //LOG_DEBUG("end parse proteoform");
}

// get mass of the modified proteoform
double Proteoform::getMass() {
  double mass = getResSeqPtr()->getSeqMass();
  for(size_t i = 0; i<change_list_.size(); i++) {
    // only unexpected changes need to to added
    if (change_list_[i]->getChangeTypePtr() == ChangeType::UNEXPECTED) {
      mass += change_list_[i]->getMassShift();
    }
  }
  return mass;
}

AlignTypePtr Proteoform::getAlignType() {
  int trunc_len = prot_mod_ptr_->getTruncPtr()->getTruncLen();
  //LOG_DEBUG("seq " << getProteinMatchSeq() << " trunc len " 
  //<< trunc_len << " start pos " << start_pos_);
  bool is_prefix = false;
  if (start_pos_ == trunc_len) {
    is_prefix = true;
  }

  bool is_suffix = false;
  if (end_pos_ == fasta_seq_ptr_->getLen() - 1) {
    is_suffix = true;
  }

  if (is_prefix) {
    if (is_suffix) {
      return AlignType::COMPLETE;
    } else {
      return AlignType::PREFIX;
    }
  } else {
    if (is_suffix) {
      return AlignType::SUFFIX;
    } else {
      return AlignType::INTERNAL;
    }
  }
}

int Proteoform::getChangeNum(ChangeTypePtr ct_ptr) {
  int n = 0;
  for (size_t i = 0; i < change_list_.size(); i++) {
    if (change_list_[i]->getChangeTypePtr() == ct_ptr) {
      n++;
    }
  }
  return n;
}

ChangePtrVec Proteoform::getChangePtrVec(ChangeTypePtr ct_ptr) {
  ChangePtrVec change_ptr_vec;
  for (size_t i = 0; i < change_list_.size(); i++) {
    if (change_list_[i]->getChangeTypePtr() == ct_ptr) {
      change_ptr_vec.push_back(change_list_[i]);
    }
  }
  return change_ptr_vec;
}

void Proteoform::addChangePtrVec(ChangePtrVec &new_changes) {
  change_list_.insert(change_list_.end(), new_changes.begin(), new_changes.end());
}

void Proteoform::addChangePtr(ChangePtr &change_ptr) {
  change_list_.push_back(change_ptr);
}

void Proteoform::rmChangePtr(ChangePtr &change_ptr) {
  for (auto iter = change_list_.begin(); iter != change_list_.end(); iter++) {
    if (*iter == change_ptr) {
      change_list_.erase(iter);
      return;
    }
  }
}

// get several segments without unexpected PTMs from a proteoform 
SegmentPtrVec Proteoform::getSegmentPtrVec() {
  ChangePtrVec changes;
  double mass_shift_sum = 0;
  for (size_t i = 0; i < change_list_.size(); i++) {
    ChangeTypePtr change_type_ptr = change_list_[i]->getChangeTypePtr();
    if (change_type_ptr == ChangeType::UNEXPECTED
        || change_type_ptr == ChangeType::VARIABLE) {
      changes.push_back(change_list_[i]);
      mass_shift_sum += change_list_[i]->getMassShift();
    }
  }
  SegmentPtrVec segments;
  double n_shift = 0;
  double c_shift = mass_shift_sum;
  int left = 0;
  for (size_t i = 0; i < changes.size(); i++) {
    int right = changes[i]->getLeftBpPos();
    SegmentPtr segment_ptr = SegmentPtr(
        new Segment(left, right, n_shift, c_shift));
    segments.push_back(segment_ptr);
    left = changes[i]->getRightBpPos();
    n_shift = n_shift + changes[i]->getMassShift();
    c_shift = c_shift - changes[i]->getMassShift();
  }
  int right = residue_seq_ptr_->getLen();
  SegmentPtr segment_ptr = SegmentPtr(
      new Segment(left, right, n_shift, c_shift));
  segments.push_back(segment_ptr);
  return segments;
}

inline void updateMatchSeq(const ChangePtrVec &changes,
                           std::vector<std::string> &left_strings,
                           std::vector<std::string> &right_strings) {
  for (size_t i = 0; i < changes.size(); i++) {
    int left_pos = changes[i]->getLeftBpPos();
    left_strings[left_pos] = "(" + left_strings[left_pos];
    int right_pos = changes[i]->getRightBpPos();
    double shift = changes[i]->getMassShift();
    right_strings[right_pos] +=  ")";
    if (ModBase::isNoneModPtr(changes[i]->getModPtr())) {
      right_strings[right_pos] = right_strings[right_pos] 
          + "["+StringUtil::convertToString(shift,5)+"]";
    } else {
      right_strings[right_pos] = right_strings[right_pos] 
          + "["+changes[i]->getModPtr()->getModResiduePtr()->getPtmPtr()->getAbbrName()+"]";
    }
  }
}
                                       
std::string Proteoform::getProteinMatchSeq() {
  std::string protein_string = fasta_seq_ptr_->getSeq();
  //LOG_DEBUG("protein string lenth " << protein_string.length() << " string " << protein_string);
  std::string mid_string = residue_seq_ptr_->toAcidString();
  //LOG_DEBUG("mid string lenth " << mid_string.length() << " string " << mid_string);
  std::sort(change_list_.begin(),change_list_.end(),Change::cmpPosInc);

  std::vector<std::string> left_strings(mid_string.length() + 1, "");
  std::vector<std::string> right_strings(mid_string.length() + 1, "");

  ChangePtrVec input_changes = getChangePtrVec(ChangeType::INPUT);
  updateMatchSeq(input_changes, left_strings, right_strings);
  ChangePtrVec fixed_changes = getChangePtrVec(ChangeType::FIXED);
  updateMatchSeq(fixed_changes, left_strings, right_strings);
  ChangePtrVec protein_var_changes = getChangePtrVec(ChangeType::PROTEIN_VARIABLE);
  updateMatchSeq(protein_var_changes, left_strings, right_strings);
  ChangePtrVec var_changes = getChangePtrVec(ChangeType::VARIABLE);
  updateMatchSeq(var_changes, left_strings, right_strings);
  ChangePtrVec unexpected_changes = getChangePtrVec(ChangeType::UNEXPECTED);
  updateMatchSeq(unexpected_changes, left_strings, right_strings);

  std::string result="";
  for (size_t i = 0; i < mid_string.length(); i++) {
    result = result + right_strings[i] + left_strings[i] + mid_string.substr(i, 1);
  }
  // last break;
  result = result + right_strings[mid_string.length()];

  std::string prefix = "";
  if(start_pos_>0) {
    prefix = protein_string.substr(start_pos_-1,1);
  }
  std::string suffix = "";
  if(end_pos_< (int)protein_string.length()-1) {
    suffix = protein_string.substr(end_pos_+1,1);
  }

  //LOG_DEBUG("Prefix " << prefix << " result " << result << " suffix length " << suffix.length() << " suffix " << suffix);
  return prefix+"."+result+"."+suffix;
}


std::string Proteoform::toString() {
  std::stringstream s;
  s<< "Begin pos: " << start_pos_ << std::endl;
  s<< "End pos: " << end_pos_ << std::endl;
  s<< "String: " << residue_seq_ptr_->toString();
  return s.str();
}

void Proteoform::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent) {
  std::string element_name = getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  fasta_seq_ptr_->appendNameDescToXml(xml_doc,element);
  prot_mod_ptr_->appendNameToXml(xml_doc,element);
  std::string str = StringUtil::convertToString(start_pos_);
  xml_doc->addElement(element, "start_pos", str.c_str());
  str = StringUtil::convertToString(end_pos_);
  xml_doc->addElement(element, "end_pos", str.c_str());
  str = StringUtil::convertToString(species_id_);
  xml_doc->addElement(element, "species_id", str.c_str());
  str = StringUtil::convertToString(prot_id_);
  xml_doc->addElement(element, "prot_id", str.c_str());

  element_name = Change::getXmlElementName() + "_list";
  xercesc::DOMElement* cl = xml_doc->createElement(element_name.c_str());
  for(size_t i=0; i<change_list_.size(); i++) {
    change_list_[i]->appendXml(xml_doc,cl);
  }
  element->appendChild(cl);
  parent->appendChild(element);
}


} /* namespace prot */

