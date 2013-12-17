/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <base/logger.hpp>

#include "base/acid.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"


namespace prot {

Acid::Acid (std::string const &name, std::string const &one_letter, 
            std::string const &three_letter, std::string const &composition, 
            double mono_mass, double avg_mass) {
  name_ = name;
  one_letter_ = one_letter;
  three_letter_ = three_letter;
  composition_ = composition;
  mono_mass_ = mono_mass;
  avg_mass_ = avg_mass;
}

AcidPtrVec getAcidPtrVecInstance(const char* file_name) {
  AcidPtrVec acid_list;
  XmlDOMParser* parser = getXmlDOMInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name);
    LOG_DEBUG( "doc " << doc);
    if (doc) {
      xercesc::DOMElement* parent = doc->getDocumentElement();
      int acid_num = getChildCount(parent, "amino_acid");
      LOG_DEBUG( "acid num " << acid_num);
      for (int i = 0; i < acid_num; i++) {
        xercesc::DOMElement* element = getChildElement(parent, "amino_acid", i);
        std::string name = getChildValue(element, "name", 0);
        std::string one_letter = getChildValue(element, "one_letter", 0);
        std::string three_letter = getChildValue(element, "three_letter", 0);
        std::string composition = getChildValue(element, "composition", 0);
        double mono_mass = getDoubleChildValue(element, "mono_mass", 0);
        double avg_mass = getDoubleChildValue(element, "average_mass", 0);
        acid_list.push_back(AcidPtr(
                new Acid(name, one_letter, three_letter, composition, mono_mass, avg_mass)));

      }
      delete doc;
    }
    // deleting parser is not necessary
  }
  return acid_list;
}

/**
 * Returns an amino acid based on the the name. Returns null if the amino
 * acid name does not exist.
 */
AcidPtr getAcidPtrByName(AcidPtrVec &acid_list, 
                         const std::string &name) {
  for (unsigned int i = 0; i < acid_list.size(); i++) {
    std::string n = acid_list[i]->getName();
    if (n == name) {
      return acid_list[i];
    }
  }
  return AcidPtr(nullptr);
}

/**
 * Returns an amino acid based on the one letter representation. Returns
 * null if the one letter representation does not exist.
 */
AcidPtr getAcidPtrByOneLetter(AcidPtrVec &acid_list, 
                              const std::string &one_letter) {
  for (unsigned int i = 0; i < acid_list.size(); i++) {
    std::string l = acid_list[i]->getOneLetter();
    if (l == one_letter)  {
      return acid_list[i];
    }
  }
  LOG_DEBUG( "Acid not found " + one_letter);
  return AcidPtr(nullptr);
}

/**
 * Returns an amino acid based on the three letter representation. Returns
 * null if the three letter representation does not exist.
 */
AcidPtr getAcidPtrByThreeLetter(AcidPtrVec &acid_list, 
                                const std::string &three_letter) {
  for (unsigned int i = 0; i < acid_list.size(); i++) {
    std::string l = acid_list[i]->getThreeLetter();
    if (l == three_letter) {
      return acid_list[i];
    }
  }
  LOG_DEBUG( "Acid not found " + three_letter);
  return AcidPtr(nullptr);
}

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool containsName(AcidPtrVec &acid_list, 
                  const std::string &name) {
  if (getAcidPtrByName(acid_list, name).get() == nullptr) {
    return false;
  }
  else {
    return true;
  }
}

/**
 * Checks if the list contains an amino acid with the specific one letter
 * representation.
 */
bool containsOneLetter(AcidPtrVec &acid_list, 
                       const std::string &one_letter) {
  if (getAcidPtrByOneLetter(acid_list, one_letter).get() == nullptr) {
    return false;
  }
  else {
    return true;
  }
}

/**
 * Checks if the list contains an amino acid with the specific three letter
 * representation.
 */
bool containsThreeLetter(AcidPtrVec &acid_vec, 
                         const std::string &three_letter) {
  if (getAcidPtrByThreeLetter(acid_vec, three_letter).get() == nullptr) {
    return false;
  }
  else {
    return true;
  }
}

/**
 * Converts a protein sequence (with one letter representation of amino
 * acids) to an amino acid array.
 */
AcidPtrVec convertSeqToAcidSeq(AcidPtrVec &acid_list, 
                               const std::string &seq) {
  AcidPtrVec acid_seq;
  if (seq.length() == 0) {
    return acid_seq;
  } else {
    for (unsigned int i = 0; i < seq.length(); i++) {
      acid_seq.push_back(
          getAcidPtrByOneLetter(acid_list, seq.substr(i, 1)));
    }
    return acid_seq;
  }
}

} /* end namespace */

