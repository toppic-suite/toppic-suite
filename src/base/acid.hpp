/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_ACID_HPP_
#define PROT_BASE_ACID_HPP_

#include <memory>
#include <string>
#include <vector>

#include "base/xml_dom_document.hpp"

namespace prot {

class Acid {
 public:
  Acid(const std::string &name, const std::string &one_letter, 
       const std::string &three_letter, const std::string &composition, 
       double mono_mass, double avg_mass); 

  Acid(xercesc::DOMElement* element); 

  void appendNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  /* Get amino acid composition. */
  std::string getAcidComposition() {return composition_;}

  /* Get average mass. */
  double getAvgMass() {return average_mass_;}

  /* Get  monoisotopic mass. */
  double getMonoMass() {return mono_mass_;}

  /* Get amino acid name. */
  std::string getName() {return name_;}

  /* Get amino acid one letter representation. */
  std::string getOneLetter() {return one_letter_;}

  /* Get amino acid three letter representation. */
  std::string getThreeLetter() {return three_letter_;}

  static std::string getXmlElementName() {return "amino_acid";}

 private:
  /* Name of amino acid */
  std::string name_;
  /* One letter representation */
  std::string one_letter_;
  /* Three letter representation */
  std::string three_letter_;
  /* amino acid chemical composition */
  std::string composition_;
  /* residue monoisotopic mass */
  double mono_mass_;
  /* residue average mass */
  double average_mass_;
};

typedef std::shared_ptr<Acid> AcidPtr;
typedef std::vector<AcidPtr> AcidPtrVec;

}

#endif

