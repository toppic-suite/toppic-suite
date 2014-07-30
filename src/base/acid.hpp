/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_ACID_HPP_
#define PROT_ACID_HPP_

#include <memory>
#include <string>
#include <vector>
#include "base/xml_dom_document.hpp"

namespace prot {

class Acid {
 public:
  Acid (const std::string &name, const std::string &one_letter, 
        const std::string &three_letter, const std::string &composition, 
        double mono_mass, double avg_mass); 

  /* Get amino acid composition. */
  std::string getAcidComposition() {return composition_;}

  /* Get average mass. */
  double getAvgMass() {return avg_mass_;}

  /* Get  monoisotopic mass. */
  double getMonoMass() {return mono_mass_;}

  /* Get amino acid name. */
  std::string getName() {return name_;}

  /* Get amino acid one letter representation. */
  std::string getOneLetter() {return one_letter_;}

  /* Get amino acid three letter representation. */
  std::string getThreeLetter() {return three_letter_;}

  void appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

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
  double avg_mass_;
};

typedef std::shared_ptr<Acid> AcidPtr;
typedef std::vector<AcidPtr> AcidPtrVec;

/* acid factory */
class AcidFactory {
 public:
  static void initFactory(const std::string &file_name);

  static const AcidPtrVec& getBaseAcidPtrVec() {return acid_ptr_vec_;}

  /**
   * Returns an amino acid based on the the name. Returns null if the amino
   * acid name does not exist.
   */
  static AcidPtr getBaseAcidPtrByName(const std::string &name);

  /**
   * Returns an amino acid based on the one letter representation. Returns
   * null if the one letter representation does not exist.
   */
  static AcidPtr getBaseAcidPtrByOneLetter(const std::string &one_letter);

  /**
   * Returns an amino acid based on the three letter representation. Returns
   * null if the three letter representation does not exist.
   */
  static AcidPtr getBaseAcidPtrByThreeLetter(const std::string &three_letter);

  /**
   * Checks if the list contains an amino acid with the specific name.
   */
  static bool baseContainsName(const std::string &name);

  /**
   * Checks if the list contains an amino acid with the specific one letter
   * representation.
   */
  static bool baseContainsOneLetter(const std::string &one_letter);

  /**
   * Checks if the list contains an amino acid with the specific three letter
   * representation.
   */
  static bool baseContainsThreeLetter(const std::string &three_letter);

  /**
   * Converts a protein sequence (with one letter representation of amino
   * acids) to an amino acid array.
   */
  static AcidPtrVec convertSeqToAcidSeq(const std::string &seq);

 private:
  static AcidPtrVec acid_ptr_vec_;
};

}
#endif

