/*
 * ptm_processor.cpp
 *
 *  Created on: Dec 20, 2013
 *      Author: xunlikun
 */

#include "ptm_processor.hpp"
#include "prsm/simple_prsm.hpp"
#include "base/fasta_reader.hpp"
#include "spec/deconv_ms.hpp"
#include "prsm/prsm.hpp"
#include "spec/spectrum_set.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "base/prot_mod.hpp"

namespace prot {

PtmProcessor::PtmProcessor(PtmMngPtr mng){
	mng_ = mng;
	init();
}

void PtmProcessor::init(){
	seqs_ = prot::readFastaToProteoform(mng_->search_db_file_name_,
                                      AcidFactory::getBaseAcidPtrVec(),
                                      ResidueFactory::getBaseResiduePtrVec(),
                                      mng_->base_data_->getDefaultProtModPtr());
	std::string sp_file_name = mng_->spectrum_file_name_;
	std::string simplePrsmFileName = mng_->spectrum_file_name_ + ".FILTER" + mng_->input_file_ext_;
	simplePrsms_  = prot::readSimplePrSM(simplePrsmFileName.c_str());
	//todo::
	prsmFindSeq(simplePrsms_,seqs_);
}

void PtmProcessor::prsmFindSeq(SimplePrSMPtrVec simple_prsms,ProteoformPtrVec seqs){
	for(unsigned int i =0;i<simple_prsms.size();i++){
		simple_prsms[i]->findSeq(seqs);
	}
}

void PtmProcessor::process(){
	PtmSearcherPtr searcher = PtmSearcherPtr(new PtmSearcher(mng_));
	processDatabase(searcher);
}

void PtmProcessor::processDatabase(PtmSearcherPtr searcher){
	std::string sp_file_name = mng_->spectrum_file_name_;
	std::string output_file_name = sp_file_name+"."+mng_->output_file_ext_;

	int n_spectra = prot::countSpNum(sp_file_name.c_str());

	MsAlignReader spReader(sp_file_name.c_str());

//	std::string output_file_name = sp_file_name+"."+mng_->output_file_ext_;
	PrSMWriterPtr all_writer= PrSMWriterPtr(new PrSMWriter(output_file_name));

	std::vector<std::vector<PrSMWriterPtr>> writers;
	for(int i=0;i<mng_->n_unknown_shift_;i++){
		std::vector<PrSMWriterPtr> temp;
		for(int j=0;j<4;j++){
			std::string file_name = output_file_name+"_"+prot::convertToString(i)+"_"+convertSemiAlignmentTypeToString(j);
			temp.push_back(PrSMWriterPtr(new PrSMWriter(file_name)));
		}
		writers.push_back(temp);
	}

	DeconvMsPtr deconv_sp;
	PrSMPtrVec3D prsms;
	for(int i=0;i<mng_->n_unknown_shift_;i++){
		PrSMPtrVec2D temp_2d;
		for(int j=0;j<4;j++){
			PrSMPtrVec temp_vec;
			for(int k=0;k<mng_->n_report_;k++){
				temp_vec.push_back(nullptr);
			}
			temp_2d.push_back(temp_vec);
		}
		prsms.push_back(temp_2d);
	}
	int cnt = 0;
	while((deconv_sp = spReader.getNextMs())!= nullptr){
		cnt++;
//		for(int i=0;i<deconv_sp->size();i++){
			double shift = prot::getProtModAcetylationShift(mng_->base_data_->getProtModPtrVec());
			SpectrumSetPtr spectrumset = prot::getSpectrumSet(deconv_sp,0,mng_->sp_para_,shift,IonTypeFactory::getBaseIonTypePtrVec());
			if(spectrumset != nullptr){
//				std::string scan = deconv_sp->getHeaderPtr()->getScansString();
				//update message;
				SimplePrSMPtrVec slectedPrsms = prot::findSimplePrsms(simplePrsms_,deconv_sp->getHeaderPtr());
				searcher->search(spectrumset,slectedPrsms,prsms);
				all_writer->writeVector3D(prsms);
				for(int j=0;j<mng_->n_unknown_shift_;j++){
					for(int k=0;k<4;k++){
						writers[j][k]->writeVector(prsms[j][k]);
					}
				}
			}
//		}
	}
	spReader.close();
}

} /* namespace prot */
