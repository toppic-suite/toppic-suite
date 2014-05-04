/*
 * xml_generator.cpp
 *
 *  Created on: Feb 24, 2014
 *      Author: xunlikun
 */

#include <map>

#include "base/file_util.hpp"
#include "prsm/prsm_species.hpp"
#include "xpp/xml_generator.hpp"

namespace prot {
XmlGenerator::XmlGenerator(std::map<std::string,std::string> arguments,
                           std::string input_file) {
  input_file_ = input_file;
  mng_ = ViewMngPtr(new ViewMng(arguments));
  anno_view_ = AnnoViewPtr(new AnnoView());
}

void XmlGenerator::outputPrsms(PrSMPtrVec prsms){
  for(unsigned int i=0;i<prsms.size();i++){
    std::string file_name = mng_->xml_path_+ FILE_SEPARATOR + 
        "prsms" + FILE_SEPARATOR + "prsm"+convertToString(prsms[i]->getId())+".xml";
    XmlWriter writer(file_name,"");
    writer.write(genePrSMView(writer.getDoc(),prsms[i], mng_->sp_para_ptr_->getMinMass()));
    writer.close();

    std::vector<std::string> file_info;
    file_info.push_back(file_name);
    file_info.push_back(mng_->executive_dir_ + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR + "prsm.xsl");
    file_info.push_back(mng_->html_path_+ FILE_SEPARATOR + "prsms" + FILE_SEPARATOR 
                        + "prsm"+ convertToString(prsms[i]->getId())+".html");
    anno_view_->file_list_.push_back(file_info);

  }

}
void XmlGenerator::outputAllPrsms(PrSMPtrVec prsms){
  std::string file_name = mng_->xml_path_+ FILE_SEPARATOR + "prsms.xml";
  XmlWriter writer(file_name,"prsm_list");
  for(unsigned int i=0;i<prsms.size();i++){
    writer.write(genePrSMView(writer.getDoc(),prsms[i], mng_->sp_para_ptr_->getMinMass()));
    writer.close();
  }
}

void XmlGenerator::outputProteins(PrSMPtrVec prsms){

  for(unsigned int i=0;i<seq_.size();i++){
    std::vector<int> species = getSpeciesIds(prsms,seq_[i]->getDbResSeqPtr()->getId());
    if(species.size()>0){
      std::string file_name = mng_->xml_path_ + FILE_SEPARATOR +"proteins" 
          +FILE_SEPARATOR+ "protein"+convertToString(seq_[i]->getDbResSeqPtr()->getId())+".xml";
      XmlWriter writer(file_name,"");
      writer.write(proteinToXml(writer.getDoc(),prsms,seq_[i],species, mng_->sp_para_ptr_->getMinMass()));
      writer.close();
      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_->executive_dir_ + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR + "protein.xsl");
      file_info.push_back(mng_->html_path_+ FILE_SEPARATOR + "proteins" + FILE_SEPARATOR 
                          + "protein"+convertToString(seq_[i]->getDbResSeqPtr()->getId())+".html");
      anno_view_->file_list_.push_back(file_info);
    }
  }
}
void XmlGenerator::outputAllProteins(PrSMPtrVec prsms){

  std::string file_name = mng_->xml_path_+ FILE_SEPARATOR +"proteins.xml";
  XmlWriter writer(file_name,"protein_list");
  writer.write(allProteinToXml(writer.getDoc(),prsms,seq_, mng_->sp_para_ptr_->getMinMass()));
  writer.close();
  std::vector<std::string> file_info;
  file_info.push_back(file_name);
  file_info.push_back(mng_->executive_dir_ + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR +"proteins.xsl");
  file_info.push_back(mng_->html_path_+ FILE_SEPARATOR + "proteins.html");
  anno_view_->file_list_.push_back(file_info);
}

void XmlGenerator::processPrSMs(PrSMPtrVec & prsms,ProteoformPtrVec proteoforms){
  MsAlignReader reader(mng_->spectrum_file_name_);
  DeconvMsPtr deconv_sp;
  int cnt = 0;
  while((deconv_sp = reader.getNextMs()) != nullptr){
    cnt++;
    for(unsigned int i=0;i<prsms.size();i++){
      if(isMatch(prsms[i],deconv_sp->getHeaderPtr())){

        //process ,deconv_sp,bpsepc
        //findseq :check proteoform's id and name

//        ProteoformPtr temp = prsms[i]->getProteoformPtr();
//        ProteoformPtr new_seq = prot::getSubProteoform(proteoforms[temp->getDbResSeqPtr()->getId()],
//                                                       temp->getStartPos(),temp->getEndPos());
//        ChangePtrVec unexpect_change_list = temp->getUnexpectedChangePtrVec();
//        new_seq->addUnexpectedChangePtrVec(unexpect_change_list);
//        prsms[i]->setProteoformPtr(new_seq);

        //initsegment:proteoform getSegmentPtrVec
        //check prsm's spectrum_scan and oriprecmass
        //refine_ms_three_ ;have read it form prsm.xml
        //initscore:have read it form prsm.xml

        prsms[i]->setDeconvMsPtr(deconv_sp);
        double delta = prsms[i]->getAdjustedPrecMass() - prsms[i]->getOriPrecMass();
        prsms[i]->setRefineMs(getMsThree(deconv_sp, delta, mng_->sp_para_ptr_));
      }
    }
  }
}
void XmlGenerator::outputSpecies(PrSMPtrVec prsms){

  std::vector<int> species = getSpeciesIds(prsms);
  for(unsigned int i=0;i<species.size();i++){
    PrSMPtrVec select_prsms = selectSpeciesPrsms(prsms,species[i]);
    if(select_prsms.size()>0){
      std::string file_name = mng_->xml_path_+ FILE_SEPARATOR + "species" 
          + FILE_SEPARATOR + "species"+convertToString(species[i])+".xml";
      XmlWriter writer(file_name,"");
      std::sort(select_prsms.begin(),select_prsms.end(),prsmEValueUp);
      writer.write(speciesToXml(writer.getDoc(),select_prsms, mng_->sp_para_ptr_->getMinMass()));
      writer.close();

      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_->executive_dir_ + FILE_SEPARATOR + "xsl" + FILE_SEPARATOR + "species.xsl");
      file_info.push_back(mng_->html_path_+ FILE_SEPARATOR+ "species" + FILE_SEPARATOR 
                          + "species"+convertToString(species[i])+".html");
      anno_view_->file_list_.push_back(file_info);
    }
  }
}
void XmlGenerator::outputFileList(){
  std::string file_name = mng_->xml_path_+ FILE_SEPARATOR + "files.xml";
  XmlWriter writer(file_name,"");
  writer.write(anno_view_->geneFileList(writer.getDoc()));
  writer.close();
}
void XmlGenerator::process(){
  prot::createFolder(mng_->xml_path_ + FILE_SEPARATOR + "species");
  prot::createFolder(mng_->xml_path_ + FILE_SEPARATOR + "prsms");
  prot::createFolder(mng_->xml_path_ + FILE_SEPARATOR + "proteins");

  ProteoformPtrVec raw_forms
        = readFastaToProteoform(mng_->database_file_name_, mng_->fix_mod_residue_list_);
  seq_ = raw_forms;
  std::string input_name = basename(mng_->spectrum_file_name_)+"."+input_file_;
  PrSMPtrVec prsms = readPrsm(basename(mng_->spectrum_file_name_)+"."+input_file_,raw_forms);
//  std::cout<<prsms[0]->getProteoformPtr()->getResSeqPtr()->toString()<<std::endl;

  processPrSMs(prsms,raw_forms);
  setSpeciesId(prsms,mng_->ppo_);
  outputPrsms(prsms);
  outputAllPrsms(prsms);
  outputSpecies(prsms);
  outputProteins(prsms);
  outputAllProteins(prsms);
  outputFileList();

}

} /* namespace prot */
