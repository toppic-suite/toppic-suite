#include "acid_vec.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace proteomics {
      
static std::vector<Acid> complete_acid_vec_;

const std::vector<Acid>* AcidVec::getInstance() {
    if (complete_acid_vec_.size() == 0) {
        proteomics::XmlDOMParser* parser = proteomics::getXmlDOMInstance();
        if (parser) {
            proteomics::XmlDOMDocument* doc = new proteomics::XmlDOMDocument(parser, "./acid.xml");
            if (doc) {
                int acid_num = doc->getChildCount("amino_acid_list", 0, "amino_acid");
                for (int i = 0; i < acid_num; i++) {
                    xercesc::DOMElement* element = doc->getElement("amino_acid", i);
                    complete_acid_vec_.push_back(proteomics::Acid(element));

                }
                delete doc;
            }
            delete parser;
        }
    }
    return &complete_acid_vec_;
}

/**
 * Returns an amino acid based on the the name. Returns null if the amino
 * acid name does not exist.
 */
Acid* AcidVec::getAcidByName(std::vector<Acid> &acid_vec, const std::string &name) {
    for (unsigned int i = 0; i < acid_vec.size(); i++) {
        std::string n = acid_vec[i].getName();
        if (n == name) {
            return &acid_vec[i];
        }
    }
    return nullptr;
}

/**
 * Returns an amino acid based on the one letter representation. Returns
 * null if the one letter representation does not exist.
 */
Acid* AcidVec::getAcidByOneLetter(std::vector<Acid> &acid_vec, const std::string &one_letter) {
    for (unsigned int i = 0; i < acid_vec.size(); i++) {
        std::string l = acid_vec[i].getOneLetter();
        if (l == one_letter) {
            return &acid_vec[i];
        }
    }
    //logger.debug("Acid not found " + one_letter);
    return nullptr;
}

/**
 * Returns an amino acid based on the three letter representation. Returns
 * null if the three letter representation does not exist.
 */
Acid* AcidVec::getAcidByThreeLetter(std::vector<Acid> &acid_vec, const std::string &three_letter) {
    for (unsigned int i = 0; i < acid_vec.size(); i++) {
        std::string l = acid_vec[i].getThreeLetter();
        if (l == three_letter) {
            return &acid_vec[i];
        }
    }
    //logger.debug("Acid not found " + three_letter);
    return nullptr;
}

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool AcidVec::containsName(std::vector<Acid> &acid_vec, const std::string &name) {
    if (getAcidByName(acid_vec, name) == nullptr) {
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
bool AcidVec::containsOneLetter(std::vector<Acid> &acid_vec, const std::string &one_letter) {
    if (getAcidByOneLetter(acid_vec, one_letter) == nullptr) {
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
bool AcidVec::containsThreeLetter(std::vector<Acid> &acid_vec, const std::string &three_letter) {
    if (getAcidByThreeLetter(acid_vec, three_letter) == nullptr) {
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
std::vector<Acid*> AcidVec::convert(std::vector<Acid> &acid_vec, const std::string &seq) {
    std::vector<Acid*> acids;
    if (seq.length() == 0) {
        return acids;
    } else {
        for (unsigned int i = 0; i < seq.length(); i++) {
            acids.push_back(AcidVec::getAcidByOneLetter(acid_vec, seq.substr(i, 1)));
        }
        return acids;
    }
}
	
};


