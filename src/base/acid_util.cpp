
#include <iostream>

#include "acid_util.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {


AcidPtrVec getAcidPtrVecInstance(const char* file_name) {
  AcidPtrVec acid_ptr_vec;
  prot::XmlDOMParser* parser = prot::getXmlDOMInstance();
  if (parser) {
    prot::XmlDOMDocument* doc = new prot::XmlDOMDocument(parser, file_name);
    if (doc) {
      int acid_num = doc->getChildCount("amino_acid_list", 0, "amino_acid");
      for (int i = 0; i < acid_num; i++) {
        xercesc::DOMElement* element = doc->getElement("amino_acid", i);
        acid_ptr_vec.push_back(AcidPtr(new Acid(element)));

      }
      delete doc;
    }
    delete parser;
  }
  return acid_ptr_vec;
}

/**
 * Returns an amino acid based on the the name. Returns null if the amino
 * acid name does not exist.
 */
AcidPtr getAcidPtrByName(AcidPtrVec &acid_ptr_vec, 
                         const std::string &name) {
  std::string x("x");
  for (unsigned int i = 0; i < acid_ptr_vec.size(); i++) {
    std::string n = acid_ptr_vec[i]->getName();
    if (n.compare(name) == 0) {
      return acid_ptr_vec[i];
    }
  }
  return AcidPtr(nullptr);
}

/**
 * Returns an amino acid based on the one letter representation. Returns
 * null if the one letter representation does not exist.
 */
AcidPtr getAcidPtrByOneLetter(AcidPtrVec &acid_ptr_vec, 
                              const std::string &one_letter) {
  for (unsigned int i = 0; i < acid_ptr_vec.size(); i++) {
    std::string l = acid_ptr_vec[i]->getOneLetter();
    if (l.compare(one_letter) == 0) {
      return acid_ptr_vec[i];
    }
  }
  //logger.debug("Acid not found " + one_letter);
  return AcidPtr(nullptr);
}

/**
 * Returns an amino acid based on the three letter representation. Returns
 * null if the three letter representation does not exist.
 */
AcidPtr getAcidPtrByThreeLetter(AcidPtrVec &acid_ptr_vec, 
                                const std::string &three_letter) {
  for (unsigned int i = 0; i < acid_ptr_vec.size(); i++) {
    std::string l = acid_ptr_vec[i]->getThreeLetter();
    if (l.compare(three_letter) == 0) {
      return acid_ptr_vec[i];
    }
  }
  //logger.debug("Acid not found " + three_letter);
  return AcidPtr(nullptr);
}

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool containsName(AcidPtrVec &acid_ptr_vec, 
                  const std::string &name) {
  if (getAcidPtrByName(acid_ptr_vec, name).get() == nullptr) {
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
bool containsOneLetter(AcidPtrVec &acid_ptr_vec, 
                       const std::string &one_letter) {
  if (getAcidPtrByOneLetter(acid_ptr_vec, one_letter).get() == nullptr) {
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
AcidPtrVec convertSeqToAcidPtrVec(AcidPtrVec &acid_ptr_vec, 
                                  const std::string &seq) {
  AcidPtrVec acid_ptr_seq;
  if (seq.length() == 0) {
    return acid_ptr_seq;
  } else {
    for (unsigned int i = 0; i < seq.length(); i++) {
      acid_ptr_seq.push_back(
          getAcidPtrByOneLetter(acid_ptr_vec, seq.substr(i, 1)));
    }
    return acid_ptr_seq;
  }
}

	
};

int main(int argc, char** argv) {
  prot::AcidPtrVec acid_ptr_vec = prot::getAcidPtrVecInstance("./acid.xml");
  std::string letter("A");
  std::cout << "A exist: " << prot::containsOneLetter(acid_ptr_vec, letter) << "\n";
  letter = "X";
  std::cout << "X exist: " << prot::containsOneLetter(acid_ptr_vec, letter) << "\n";
  exit(0);
}


