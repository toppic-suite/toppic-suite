/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

#include "base/ptm.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

PtmPtrVec PtmFactory::ptm_ptr_vec_;
PtmPairVec PtmFactory::ptm_pair_vec_;

Ptm::Ptm(const std::string &abbr_name, double mono_mass):
    abbr_name_(abbr_name), mono_mass_(mono_mass) {
}

Ptm::Ptm(const std::string &name, const std::string &abbr_name,
         double mono_mass, const std::string &posN,
         const std::string &posC, const std::string & pos, int id):
    name_(name), abbr_name_(abbr_name), mono_mass_(mono_mass),
    posN_(posN), posC_(posC), pos_(pos),
    unimod_id_(id) {
}

void Ptm::appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent) {
    xercesc::DOMElement* element = xml_doc->createElement("modification");
    xml_doc->addElement(element, "abbr_name", abbr_name_.c_str());
    std::string str = convertToString(mono_mass_);
    xml_doc->addElement(element, "mono_mass", str.c_str());
    str = convertToString(unimod_id_);
    xml_doc->addElement(element, "unimod", str.c_str());
    parent->appendChild(element);
}

void PtmFactory::initFactory(const std::string &file_name) {
    XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
    if (parser) {
        XmlDOMDocument doc(parser, file_name.c_str());
        xercesc::DOMElement* parent = doc.getDocumentElement();
        int ptm_num = getChildCount(parent, "modification");
        for (int i = 0; i < ptm_num; i++) {
            xercesc::DOMElement* element = getChildElement(parent, "modification", i);
            std::string name = getChildValue(element, "name", 0);
            std::string abbr_name = getChildValue(element, "abbreviation", 0);
            int id = std::stoi(getChildValue(element, "unimod", 0));
            std::string posN, posC, pos;
            double mono_mass = getDoubleChildValue(element, "mono_mass", 0);
            xercesc::DOMNodeList* list = element->getElementsByTagName(X("residues"));
            for (size_t j = 0; j < list->getLength(); j++) {
                if (Y(list->item(j)->getAttributes()->item(0)->getNodeValue()).compare("N-term") == 0) {
                    xercesc::DOMElement* child = dynamic_cast<xercesc::DOMElement*>(list->item(j));
                    posN = Y(child->getTextContent());
                } else if (Y(list->item(j)->getAttributes()->item(0)->getNodeValue()).compare("C-term") == 0) {
                    xercesc::DOMElement* child = dynamic_cast<xercesc::DOMElement*>(list->item(j));
                    posC = Y(child->getTextContent());
                } else if (Y(list->item(j)->getAttributes()->item(0)->getNodeValue()).compare("Anywhere") == 0) {
                    xercesc::DOMElement* child = dynamic_cast<xercesc::DOMElement*>(list->item(j));
                    pos = Y(child->getTextContent());
                }
            }
            ptm_ptr_vec_.push_back(
                PtmPtr(new Ptm(name, abbr_name, mono_mass, posN, posC, pos, id)));

        }
    }

    std::sort(ptm_ptr_vec_.begin(), ptm_ptr_vec_.end(), ptmMassUp);

    for (size_t i = 1; i < ptm_ptr_vec_.size(); i++) {
        if (ptm_ptr_vec_[i]->getName() == "Carbamidomethylation" 
            || ptm_ptr_vec_[i]->getName() == "Carboxymethyl"
                || ptm_ptr_vec_[i]->getMonoMass() == 0.0)
            continue;
        for (size_t j = 1; j < ptm_ptr_vec_.size(); j++) {
            if (ptm_ptr_vec_[j]->getName() == "Carbamidomethylation" 
                || ptm_ptr_vec_[j]->getName() == "Carboxymethyl"
                    || ptm_ptr_vec_[j]->getMonoMass() == 0.0)
                continue;

            ptm_pair_vec_.push_back(std::make_pair(ptm_ptr_vec_[i], ptm_ptr_vec_[j]));
        }
    }

    std::sort(ptm_pair_vec_.begin(), ptm_pair_vec_.end(), ptmPairMassUp);

}

PtmPtr PtmFactory::findEmptyPtmPtr() {
    for (size_t i = 0; i < ptm_ptr_vec_.size(); i++) {
        if (ptm_ptr_vec_[i]->isEmpty()) {
            return ptm_ptr_vec_[i];
        }
    }
    throw "Empty ptm does not exist!";
}

/**
 *   Returns a PTM based on the abbreviation name. Returns null if the
 *   abbreviation name does not exist.
 */
PtmPtr PtmFactory::getBasePtmPtrByAbbrName(const std::string &abbr_name) {
    for (size_t i = 0; i < ptm_ptr_vec_.size(); i++) {
        std::string n = ptm_ptr_vec_[i]->getAbbrName();
        if (n == abbr_name) {
            return ptm_ptr_vec_[i];
        }
    }
    return PtmPtr(nullptr);
}

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool PtmFactory::baseContainAbbrName(const std::string &abbr_name) {
    return getBasePtmPtrByAbbrName(abbr_name).get() != nullptr;
}

PtmPtr PtmFactory::addBasePtm(const std::string &abbr_name, double mono_mass) {
    PtmPtr ptm_ptr = getBasePtmPtrByAbbrName(abbr_name);
    if (ptm_ptr.get() == nullptr) {
        PtmPtr new_ptm(new Ptm(abbr_name, mono_mass));
        ptm_ptr_vec_.push_back(new_ptm);
        return new_ptm;
    } else {
        return ptm_ptr;
    }
}

PtmPtrVec PtmFactory::getBasePtmPtrByMass(double mass, double error_tolerance) {
    PtmPtrVec res;
    for (size_t i = 0; i < ptm_ptr_vec_.size(); i++) {

        if (ptm_ptr_vec_[i]->getMonoMass() == 0.0) {
            continue;
        }

        double ptm_mass = ptm_ptr_vec_[i]->getMonoMass();
        if (std::abs(ptm_mass - mass) < error_tolerance
                || std::abs(std::abs(ptm_mass - mass) - 1) < error_tolerance)
            res.push_back(ptm_ptr_vec_[i]);

    }
    return res;
}

PtmPairVec PtmFactory::getBasePtmPairByMass(double mass1, double mass2,
        double error_tolerance) {
    PtmPairVec res;
    double mass = mass1 + mass2;
    for (size_t i = 0; i < ptm_pair_vec_.size(); i++) {
        double pair_mass = ptm_pair_vec_[i].first->getMonoMass()
                           + ptm_pair_vec_[i].second->getMonoMass();
        if (std::abs(pair_mass - mass) < error_tolerance
                || std::abs(std::abs(pair_mass - mass) - 1) < error_tolerance)
            res.push_back(ptm_pair_vec_[i]);

    }
    return res;
}

bool PtmFactory::isKnown(double m, double error_tolerance) {
    for (size_t i = 0; i < ptm_ptr_vec_.size(); i++) {
        if (ptm_ptr_vec_[i]->getMonoMass() == 0.0
            ||ptm_ptr_vec_[i]->getAbbrName() == "Carbamidomethylation" 
                || ptm_ptr_vec_[i]->getAbbrName() == "Carboxymethyl" ) {
            continue;
        }
        double ptm_mass = ptm_ptr_vec_[i]->getMonoMass();
        if (std::abs(ptm_mass - m) < error_tolerance
                || std::abs(std::abs(ptm_mass - m) - 1) < error_tolerance) {
            return true;
        }
    }
    return false;
}

}

