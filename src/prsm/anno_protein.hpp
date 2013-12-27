/*
 * anno_protein.hpp
 *
 *  Created on: Dec 26, 2013
 *      Author: xunlikun
 */

#ifndef ANNO_PROTEIN_HPP_
#define ANNO_PROTEIN_HPP_

#include <xercesc/dom/DOM.hpp>
#include "prsm/sem_align_type.hpp"
#include "prsm/anno_cleavage.hpp"
#include "base/proteoform.hpp"
#include "base/ion_type.hpp"
#include "base/mass_constant.hpp"
#include "spec/extend_peak.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

class AnnoProtein {
public:
	AnnoProtein(ProteoformPtr seq,int firstResPos,int lastResPos,bool nAcetylation,int nKnownShifts);
	AnnoProtein(xercesc::DOMElement * element);
	double getMassShiftSum(){
		double mass = 0;
		ChangePtrVec changes = seq_->getChangePtrVec();
		for(int i=0;i<changes.size();i++){
			mass += changes[i]->getMassShift();
		}
		return mass;
	}

	int getMassShiftNum(){
		return seq_->getChangePtrVec().size();
	}

	double getMass(){
		std::vector<double> extBMass = seq_->getBpSpecPtr()->getBreakPointMasses(IonTypePtr(new IonType("B",true,0)));
		double mass = extBMass[last_res_pos_+1]-extBMass[first_res_pos_]+ MassConstant::getWaterMass();
		return mass + getMassShiftNum();
	}

	ProteoformPtr getProteoformPtr(){return seq_;}
	int getFirstResPos(){return first_res_pos_;}
	int getLastResPos(){return last_res_pos_;}
	int getAnnoSeqLen(){return last_res_pos_ - first_res_pos_;}
	bool isnAcetylation(){return n_acetylation_;}
	std::string getSeqName(){return seq_name_;}
	int getSeqId(){return seq_id_;}
	void setSpeciesId(int sepecies_id){species_id_ = sepecies_id;}
	int getSpeciesId(){return species_id_;}
	bool isNTermShift(){
		if(seq_->getSegmentPtrVec()[0]->getPepNTermShift()==0){
			return false;
		}
		return true;
	}
	void setNKnownShifts(int n_known_shifts){n_expected_shifts_=n_known_shifts;}
	int getKnowShiftNum(){return n_expected_shifts_;}
	std::string getScan(){return scan_;}

	int getUnKnowShiftNum();
	double getTRuncSeqMass();
	SemiAlignTypePtr getAlignType(TruncPtrVec allowedNTruncs, TruncPtrVec allowedCTruncs,double thresh);
	void process();
	void findSeq(ProteoformPtrVec seqs);
	std::string getAnnoMatchSeq();
	ResiduePtrVec getAnnoResidues();
	AnnoCleavagePtrVec getAnnoCleavage(ExtendMsPtr ms_three,double min_mass);
	TheoPeakPtrVec getTheoPeaks(ActivationPtr type,double min_mass);
	PeakIonPairPtrVec getMatchedPairs(PeakIonPairPtrVec pairs,int peakId);
	xercesc::DOMElement * ionToXml();
	xercesc::DOMElement * peakToXml();

private:
	int seq_id_;
	std::string seq_name_;
	int first_res_pos_;
	int last_res_pos_;
	int n_expected_shifts_;
	bool n_acetylation_;
	int species_id_;
	ProteoformPtr seq_;
	std::string scan_;

	void init();
	void initSegments();
};

} /* namespace prot */

#endif /* ANNO_PROTEIN_HPP_ */
