/*
 * basic_diag_pair.cpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

#include <ptmsearch/basic_diag_pair.hpp>

namespace prot {

BasicDiagPairDiagPtrVec getDiagonals(DiagonalHeaderPtrVec headers,PrmMsPtr ms_six,ProteoformPtr seq,PtmMngPtr mng){
	BasicDiagPairDiagPtrVec dianonal_list;
	int cnt =0;
	for(int i=0;i<headers.size();i++){
		BasicDiagPairDiagPtr diagonal =
	}
}

BasicDiagPairDiagPtr getDiagonal(int cnt,DiagonalHeaderPtr header,PrmMsPtr ms_six,ProteoformPtr seq,PtmMngPtr mng){
	double n_shift = header->getProtNTermShift();
	double c_shift = ms_six->getHeaderPtr()->getPrecMonoMass()-seq->getBpSpecPtr()->getResSeqMass()-n_shift;

}
} /* namespace prot */
