/*
 * xml_generator.cpp
 *
 *  Created on: Feb 24, 2014
 *      Author: xunlikun
 */

#include <xpp/xml_generator.hpp>

namespace prot {
XmlGenerator::XmlGenerator(std::string spec_file,std::string db_file,std::string input_file){
  spec_file_ = spec_file;
  db_file_=db_file;
  input_file_ = input_file;
  output_file_ = "xml/";
}
void XmlGenerator::outputPrsms(PrSMPtrVec prsms){
  for(unsigned int i=0;i<prsms.size();i++){
    std::string file_name = output_file_+"prsm"+convertToString(prsms[i]->getId())+".xml";
    XmlWriter writer(file_name,"");
    writer.write(genePrSMView(writer.getDoc(),prsms[i]));
  }
}
void XmlGenerator::outputAllPrsms(PrSMPtrVec prsms){
  std::string file_name = output_file_+"prsms.xml";
  XmlWriter writer(file_name,"prsm_list");
  for(unsigned int i=0;i<prsms.size();i++){
    writer.write(genePrSMView(writer.getDoc(),prsms[i]));
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
void XmlGenerator::process(){
  ProteoformPtrVec proteoforms = prot::readFastaToProteoform(db_file_,ResidueFactory::getBaseResiduePtrVec());
  PrSMPtrVec prsms = readPrsm(basename(spec_file_)+"."+input_file_,proteoforms);
  outputPrsms(prsms);
  outputAllPrsms(prsms);
  outputProteins(prsms);
  outputAllProteins(prsms);
}
} /* namespace prot */
