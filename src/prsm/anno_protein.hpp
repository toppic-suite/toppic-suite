/*
 * anno_protein.hpp
 *
 *  Created on: Dec 26, 2013
 *      Author: xunlikun
 */

#ifndef ANNO_PROTEIN_HPP_
#define ANNO_PROTEIN_HPP_

namespace prot {

class AnnoProtein {
public:
	AnnoProtein();
	int seq_id_;
	std::string seq_name_;
	int first_res_pos_;
	int last_res_pos_;
	int n_expected_shifts_;
	bool n_acetylation_;
	int species_id_;
	ProteoformPtr seq_;
	std::string scan_;

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
		return seq_->getBpSpecPtr()->getBreakPointMasses(IonTypePtr(new IonType("B",true,0)));
	}

	std::string getAnnoMatchSeq(){

	}

};

} /* namespace prot */

#endif /* ANNO_PROTEIN_HPP_ */
