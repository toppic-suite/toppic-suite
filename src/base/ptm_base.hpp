/*
 * author  Xiaowen Liu
 * date    2013-11-17
 */

#ifndef PROT_PTM_HPP_
#define PROT_PTM_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

#define PTM_ACETYLATION "Acetylation"

class Ptm;
typedef std::shared_ptr<Ptm> PtmPtr;
typedef std::vector<PtmPtr> PtmPtrVec;
typedef std::pair<PtmPtr, PtmPtr> PtmPair;
typedef std::vector<PtmPair> PtmPairVec;

class Ptm {
  public:
    Ptm(const std::string &abbr_name, double mono_mass);

    Ptm(const std::string &name, const std::string &abbr_name,
        double mono_mass, const std::string &posN,
        const std::string &posC, const std::string &pos, int id);

    const std::string& getAbbrName() {
        return abbr_name_;
    }

    const std::string& getName() {
        return name_;
    }

    /* Get  monoisotopic mass. */
    double getMonoMass() {
        return mono_mass_;
    }

    /* Is it an empty PTM. */
    bool isEmpty() {
        return mono_mass_ == 0.0;
    };

    /* Is it the PTM acetylation */
    bool isAcetylation() {
        return abbr_name_ == PTM_ACETYLATION;
    };

    std::string getPosN() {
        return posN_;
    }

    std::string getPosC() {
        return posC_;
    }

    std::string getPos() {
        return pos_;
    }

    int getID() {
        return unimod_id_;
    }

    void appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  private:
    /* Full name and abbreviation name of a PTM */
    std::string name_, abbr_name_;
    /* monoisotopic mass */
    double mono_mass_;
    // possible positions
    std::string posN_, posC_, pos_;
    // unimod id
    int unimod_id_;
};

/* ptm factory */
class PtmFactory {
  public:
    static void initFactory(const std::string &file_name);
    static const PtmPtrVec& getBasePtmPtrVec() {
        return ptm_ptr_vec_;
    }
    static void selectPtm(const std::string & file_name);
    static PtmPtr findEmptyPtmPtr();
    /**
     * Returns a PTM based on the abbreviation name. Returns null if the
     * abbreviation name does not exist.
     */
    static PtmPtr getBasePtmPtrByAbbrName(const std::string &abbr_name);
    /**
     * Checks if the list contains an amino acid with the specific name.
     */
    static bool baseContainAbbrName(const std::string &abbr_name);

    static PtmPtr addBasePtm(const std::string &abbr_name, double mono_mass);

    static PtmPtr getPtmPtr_Acetylation() {
        return getBasePtmPtrByAbbrName(PTM_ACETYLATION);
    }

    static PtmPtrVec getBasePtmPtrByMass(double mass, double error_tolerance);
    static PtmPairVec getBasePtmPairByMass(double mass1, double mass2,
                                           double error_tolerance);

    static bool isKnown(double m, double error_tolerance);

  private:
    static PtmPtrVec ptm_ptr_vec_;
    static PtmPairVec ptm_pair_vec_;
    static PtmPtrVec ptm_select_vec_;
    static std::vector<std::string> ptm_select_;
};

inline bool ptmMassUp(const PtmPtr &a, const PtmPtr &b) {
    return a->getMonoMass() < b->getMonoMass();
}

inline bool ptmPairMassUp(const PtmPair &a, const PtmPair &b) {
    return a.first->getMonoMass() + a.second->getMonoMass()
           < b.first->getMonoMass() + b.second->getMonoMass();
}

}
#endif

