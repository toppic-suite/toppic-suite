/*
 * basic_diag_pair.cpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

#include <ptmsearch/basic_diag_pair.hpp>

namespace prot {

BasicDiagPairDiagPtrVec getDiagonals(DiagonalHeaderPtrVec headers,PrmMsPtr ms_six,ProteoformPtr seq,PtmMngPtr mng){
	BasicDiagPairDiagPtrVec diagonal_list;
	int cnt =0;
	for(int i=0;i<headers.size();i++){
		BasicDiagPairDiagPtr diagonal = getDiagonal(cnt,headers[i],ms_six,seq,mng);
		if(diagonal!=nullptr){
			diagonal_list.push_back(diagonal);
			cnt++;
		}
	}
	return diagonal_list;
}

BasicDiagPairDiagPtr getDiagonal(int cnt,DiagonalHeaderPtr header,PrmMsPtr ms_six,ProteoformPtr seq,PtmMngPtr mng){
	double n_shift = header->getProtNTermShift();
	double c_shift = ms_six->getHeaderPtr()->getPrecMonoMass()-seq->getBpSpecPtr()->getResSeqMass()-n_shift;
	prot::setPrefixSuffix(header,c_shift,seq,mng);
	BasicDiagPairPtrVec diag_pair_list = prot::compDiagPair(ms_six,seq->getBpSpecPtr()->getBreakPointMasses(IonTypePtr(new IonType("B",true,0))),header);
	if(diag_pair_list.size()>0){
		header->setId(cnt);
		return BasicDiagPairDiagPtr(new Diagonal<BasicDiagPair>(diag_pair_list,header));
	}
	return nullptr;
}
} /* namespace prot */
