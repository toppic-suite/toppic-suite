/*
 * xml_generator.cpp
 *
 *  Created on: Feb 24, 2014
 *      Author: xunlikun
 */

#include <map>
#include <base/species.hpp>
#include <xpp/xml_generator.hpp>

namespace prot {
XmlGenerator::XmlGenerator(std::map<std::string,std::string> arguments,std::string input_file){
  spec_file_ = arguments["spectrumFileName"];
  db_file_=arguments["databaseFileName"];
  arguments_ = arguments;
  input_file_ = input_file;
  output_file_ = "xml/";
  ppo_ = atoi(arguments["errorTolerance"].c_str())*0.000001;
  PtmMngPtr ptm_search_mng = PtmMngPtr(new PtmMng(arguments));
  sp_para_ptr_ = ptm_search_mng->sp_para_;
  mng_ = ViewMngPtr(new ViewMng(arguments["configuration"]));
}
void XmlGenerator::outputPrsms(PrSMPtrVec prsms){
  for(unsigned int i=0;i<prsms.size();i++){
    std::string file_name = output_file_+"prsm"+convertToString(prsms[i]->getId())+".xml";
    XmlWriter writer(file_name,"");
    writer.write(prot::genePrSMView(writer.getDoc(),prsms[i]));
  }
}
void XmlGenerator::outputAllPrsms(PrSMPtrVec prsms){
  std::string file_name = output_file_+"prsms.xml";
  XmlWriter writer(file_name,"prsm_list");
  for(unsigned int i=0;i<prsms.size();i++){
    writer.write(prot::genePrSMView(writer.getDoc(),prsms[i]));
  }
}
void XmlGenerator::outputProteins(PrSMPtrVec prsms){
  for(unsigned int i=0;i<prsms.size();i++){
    std::string file_name = output_file_+"protein"+convertToString(prsms[i]
                                                                         ->getProteoformPtr()
                                                                         ->getDbResSeqPtr()
                                                                         ->getId())+".xml";
    XmlWriter writer(file_name,"");
    writer.write(geneProteinView(writer.getDoc(),prsms[i]->getProteoformPtr(),
                                 prsms[i]->getRefineMs(),prsms[i]->getMinMass()));
  }
}
void XmlGenerator::outputAllProteins(PrSMPtrVec prsms){
  std::string file_name = output_file_+"proteins.xml";
  XmlWriter writer(file_name,"protein_list");
  for(unsigned int i=0;i<prsms.size();i++){
    writer.write(geneProteinView(writer.getDoc(),prsms[i]->getProteoformPtr(),
                                 prsms[i]->getRefineMs(),prsms[i]->getMinMass()));
  }
}

//bool isMatch(PrSMPtr & prsm,MsHeaderPtr header){
//  int id = header->getId();
//  std::string scan = header->getScansString();
//  int prec_id = header->getPrecId();
//  double prec_mass = header->getPrecMonoMass();
//  if(id==prsm->getSpectrumId() && prec_id==prsm->getPrecurorId()){
//    if(scan.compare(prsm->getSpectrumScan())!=0 ||prec_mass != prsm->getOriPrecMass()){
//      std::cout<<"Error in Prsm"<<std::endl;
//    }
//    return true;
//  }
//  return false;
//}
void XmlGenerator::processPrSMs(PrSMPtrVec & prsms,ProteoformPtrVec proteoforms){
  MsAlignReader reader(spec_file_);
  DeconvMsPtr deconv_sp;
  int cnt = 0;
  while((deconv_sp = reader.getNextMs()) != nullptr){
    cnt++;
    for(unsigned int i=0;i<prsms.size();i++){
      if(isMatch(prsms[i],deconv_sp->getHeaderPtr())){

        //process ,deconv_sp,bpsepc
        //findseq :check proteoform's id and name
        ProteoformPtr temp = prsms[i]->getProteoformPtr();
        ProteoformPtr new_seq = prot::getSubProteoform(proteoforms[temp->getDbResSeqPtr()->getId()],
                                                       temp->getStartPos(),temp->getEndPos());
        ChangePtrVec unexpect_change_list = temp->getUnexpectedChangePtrVec();
        new_seq->addUnexpectedChangePtrVec(unexpect_change_list);
//        new_seq->setResSeqPtr(temp->getResSeqPtr());
        prsms[i]->setProteoformPtr(new_seq);

        //initsegment:proteoform getSegmentPtrVec
        //check prsm's spectrum_scan and oriprecmass
        //refine_ms_three_ ;have read it form prsm.xml
        //initscore:have read it form prsm.xml

        prsms[i]->setDeconvMsPtr(deconv_sp);
        double delta = prsms[i]->getAdjustedPrecMass() - prsms[i]->getOriPrecMass();
        prsms[i]->setRefineMs(getMsThree(deconv_sp, delta, sp_para_ptr_));
      }
    }
  }
}
void XmlGenerator::process(){
  ProteoformPtrVec raw_forms
        = readFastaToProteoform(db_file_,
                                mng_->base_data_ptr_->getFixModResiduePtrVec());
  ProteoformPtrVec proteoforms
        = generateProtModProteoform(raw_forms,
                                    mng_->base_data_ptr_->getAllowProtModPtrVec());
  PrSMPtrVec prsms = readPrsm(basename(spec_file_)+"."+input_file_,raw_forms);
  processPrSMs(prsms,raw_forms);
  setSpeciesId(prsms,ppo_);
  outputPrsms(prsms);
  outputAllPrsms(prsms);
  outputProteins(prsms);
  outputAllProteins(prsms);
}

} /* namespace prot */
