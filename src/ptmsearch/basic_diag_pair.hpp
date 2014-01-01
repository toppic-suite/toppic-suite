/*
 * basic_diag_pair.hpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

#ifndef BASIC_DIAG_PAIR_HPP_
#define BASIC_DIAG_PAIR_HPP_

#include <memory>
#include <vector>

namespace prot {
class BasicDiagPair;
typedef std::shared_ptr<BasicDiagPair> BasicDiagPairPtr;
typedef std::vector<BasicDiagPairPtr> BasicDiagPairPtrVec;
typedef std::shared_ptr<Diagonal<BasicDiagPairPtr>> BasicDiagPairDiagPtr;
class BasicDiagPair:public Pair {
public:
	BasicDiagPair(int x,int y,double score,int diag_order,double diff,int prm_peak_type);
	BasicDiagPair(int x,int y,double score,int diag_order,double diff,int prm_peak_type,BasicDiagPairDiagPtr diagonal);
	BasicDiagPair(BasicDiagPairPtr pair);

	int getBaseType() const {
		return base_type_;
	}

	int getDiagOrder() const {
		return diag_order_;
	}

	const BasicDiagPairDiagPtr& getDiagonal() const {
		return diagonal_;
	}

	void setDiagonal(const BasicDiagPairDiagPtr& diagonal) {
		diagonal_ = diagonal;
	}

	double getDiff() const {
		return diff_;
	}

	double getScore() const {
		return score_;
	}

protected:
	int diag_order_;
	double diff_;
	BasicDiagPairDiagPtr diagonal_;
	double score_;
	int base_type_;

};

BasicDiagPairPtrVec compDiagPair(PrmMsPtr sp,std::vector<double> seq_masses,DiagonalHeaderPtr header);
bool contains(BasicDiagPairPtrVec pairs,int y){
	for(int i=0;i<pairs.size();i++){
		if(y==pairs[i]->getY()){
			return true;
		}
	}
	return false;
}
} /* namespace prot */

#endif /* BASIC_DIAG_PAIR_HPP_ */
