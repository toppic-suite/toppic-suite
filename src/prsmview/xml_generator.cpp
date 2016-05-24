#include <map>

#include "base/file_util.hpp"
#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_species.hpp"
#include "prsmview/anno_prsm.hpp"
#include "prsmview/anno_view.hpp"
#include "prsmview/xml_generator.hpp"

namespace prot {
XmlGenerator::XmlGenerator(PrsmParaPtr prsm_para_ptr, 
                           const std::string &exec_dir, 
                           const std::string &input_file_ext) {
  input_file_ext_ = input_file_ext;
  mng_ptr_ = PrsmViewMngPtr(new PrsmViewMng(prsm_para_ptr, exec_dir));
  anno_view_ptr_ = AnnoViewPtr(new AnnoView());
}

void XmlGenerator::outputPrsms(const PrsmPtrVec &prsm_ptrs){
  for(size_t i = 0; i < prsm_ptrs.size(); i++){
    mng_ptr_->cnt_++;
    std::string file_name = mng_ptr_->xml_path_+ FileUtil::getFileSeparator() + 
        "prsms" + FileUtil::getFileSeparator() + "prsm"+StringUtil::convertToString(prsm_ptrs[i]->getPrsmId())+".xml";
    XmlWriter writer(file_name,"");
    std::cout << std::flush << "Generating xml files - processing " << mng_ptr_->cnt_ << " of " << mng_ptr_->num_files_ << " files.\r";
    writer.write(geneAnnoPrsm(writer.getDoc(),prsm_ptrs[i], mng_ptr_));
    writer.close();

    LOG_DEBUG("output prsm completed " << i );

    std::vector<std::string> file_info;
    file_info.push_back(file_name);
    file_info.push_back(mng_ptr_->executive_dir_ + FileUtil::getFileSeparator() + "toppic_resources" 
                        + FileUtil::getFileSeparator() + "xsl" + FileUtil::getFileSeparator() + "prsm.xsl");
    file_info.push_back(mng_ptr_->html_path_+ FileUtil::getFileSeparator() + "prsms" + FileUtil::getFileSeparator() 
                        + "prsm"+ StringUtil::convertToString(prsm_ptrs[i]->getPrsmId())+".html");
    anno_view_ptr_->file_list_.push_back(file_info);
  }
}

void XmlGenerator::outputAllPrsms(const PrsmPtrVec &prsm_ptrs){
  std::string file_name = mng_ptr_->xml_path_+ FileUtil::getFileSeparator() + "prsms.xml";
  XmlWriter writer(file_name,"prsm_list");
  for(unsigned int i=0;i<prsm_ptrs.size();i++){
    mng_ptr_->cnt_++;
    std::cout << std::flush << "Generating xml files - processing " << mng_ptr_->cnt_ << " of " << mng_ptr_->num_files_ << " files.\r";
    writer.write(geneAnnoPrsm(writer.getDoc(),prsm_ptrs[i], mng_ptr_));
    writer.close();
  }
  std::cout << std::endl;
}

void XmlGenerator::outputProteoforms(const PrsmPtrVec &prsm_ptrs){

  std::vector<int> species_ids = PrsmUtil::getSpeciesIds(prsm_ptrs);
  LOG_DEBUG("species id size " << species_ids.size());
  for(unsigned int i=0;i<species_ids.size();i++){
    mng_ptr_->cnt_++;
    std::cout << std::flush << "Generating xml files - processing " << mng_ptr_->cnt_ << " of " << mng_ptr_->num_files_ << " files.\r";
    PrsmPtrVec select_prsm_ptrs = PrsmUtil::selectSpeciesPrsms(prsm_ptrs,species_ids[i]);
    if(select_prsm_ptrs.size()>0){
      std::string file_name = mng_ptr_->xml_path_+ FileUtil::getFileSeparator() + "proteoforms" 
          + FileUtil::getFileSeparator() + "proteoform"+StringUtil::convertToString(species_ids[i])+".xml";
      XmlWriter writer(file_name,"");
      std::sort(select_prsm_ptrs.begin(),select_prsm_ptrs.end(),Prsm::cmpEValueInc);
      writer.write(proteoformToXml(writer.getDoc(),select_prsm_ptrs, mng_ptr_));
      writer.close();
      LOG_DEBUG("output proteoform completed " << i);

      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_ptr_->executive_dir_ + FileUtil::getFileSeparator() + 
                          "toppic_resources" + FileUtil::getFileSeparator() + "xsl" + FileUtil::getFileSeparator() + "proteoform.xsl");
      file_info.push_back(mng_ptr_->html_path_+ FileUtil::getFileSeparator()+ "proteoforms" + FileUtil::getFileSeparator() 
                          + "proteoform"+StringUtil::convertToString(species_ids[i])+".html");
      anno_view_ptr_->file_list_.push_back(file_info);
    }
  }
}

void XmlGenerator::outputProteins(const PrsmPtrVec &prsm_ptrs){
  //LOG_DEBUG("prsm number " << prsm_ptrs.size());
  FastaReader reader(mng_ptr_->prsm_para_ptr_->getSearchDbFileName());
  FastaSeqPtr seq_ptr = reader.getNextSeq();

  while (seq_ptr != nullptr) {
    mng_ptr_->cnt_++;
    std::cout << std::flush << "Generating xml files is processing " << mng_ptr_->cnt_ << " of " << mng_ptr_->num_files_ << " files.\r";
    std::string seq_name = seq_ptr->getName();
    int prot_id = PrsmUtil::getProteinId(prsm_ptrs, seq_name);
    std::vector<int> species = PrsmUtil::getSpeciesIds(prsm_ptrs,seq_name);
    //LOG_DEBUG("species size " << species.size());
    if(species.size()>0){
      std::string file_name = mng_ptr_->xml_path_ + FileUtil::getFileSeparator() +"proteins" 
          +FileUtil::getFileSeparator()+ "protein"+StringUtil::convertToString(prot_id)+".xml";
      XmlWriter writer(file_name,"");
      writer.write(proteinToXml(writer.getDoc(),prsm_ptrs,seq_ptr,prot_id, species, mng_ptr_));
      writer.close();
      LOG_DEBUG("output protein completed " << seq_name);
      std::vector<std::string> file_info;
      file_info.push_back(file_name);
      file_info.push_back(mng_ptr_->executive_dir_ + FileUtil::getFileSeparator() + "toppic_resources" 
                          + FileUtil::getFileSeparator() + "xsl" + FileUtil::getFileSeparator() + "protein.xsl");
      file_info.push_back(mng_ptr_->html_path_+ FileUtil::getFileSeparator() + "proteins" + FileUtil::getFileSeparator() 
                          + "protein"+StringUtil::convertToString(prot_id)+".html");
      anno_view_ptr_->file_list_.push_back(file_info);
    }
    seq_ptr = reader.getNextSeq();
  }
}

void XmlGenerator::outputAllProteins(const PrsmPtrVec &prsm_ptrs){

  std::string file_name = mng_ptr_->xml_path_+ FileUtil::getFileSeparator() +"proteins.xml";
  XmlWriter writer(file_name,"protein_list");
  writer.write(allProteinToXml(writer.getDoc(), prsm_ptrs, mng_ptr_));
  writer.close();
  std::vector<std::string> file_info;
  file_info.push_back(file_name);
  file_info.push_back(mng_ptr_->executive_dir_ + FileUtil::getFileSeparator() + 
                      "toppic_resources" + FileUtil::getFileSeparator() + "xsl" + FileUtil::getFileSeparator() +"proteins.xsl");
  file_info.push_back(mng_ptr_->html_path_+ FileUtil::getFileSeparator() + "proteins.html");
  anno_view_ptr_->file_list_.push_back(file_info);
}

void XmlGenerator::outputFileList(){
  std::string file_name = mng_ptr_->xml_path_+ FileUtil::getFileSeparator() + "files.xml";
  XmlWriter writer(file_name,"");
  writer.write(anno_view_ptr_->geneFileList(writer.getDoc()));
  writer.close();
}
void XmlGenerator::process(){
  LOG_DEBUG("process start");
  FileUtil::createFolder(mng_ptr_->xml_path_ + FileUtil::getFileSeparator() + "proteoforms");
  FileUtil::createFolder(mng_ptr_->xml_path_ + FileUtil::getFileSeparator() + "prsms");
  FileUtil::createFolder(mng_ptr_->xml_path_ + FileUtil::getFileSeparator() + "proteins");
  LOG_DEBUG("fold created");

  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string spectrum_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = FileUtil::basename(spectrum_file_name) + "." + input_file_ext_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();

  PrsmPtrVec prsm_ptrs = PrsmReader::readAllPrsms(input_file_name, db_file_name, fix_mod_ptr_vec);

  LOG_DEBUG("prsm loaded");

  PrsmUtil::addSpectrumPtrsToPrsms(prsm_ptrs, prsm_para_ptr);
  LOG_DEBUG("spectrum added");

  std::sort(prsm_ptrs.begin(), prsm_ptrs.end(), Prsm::cmpEValueInc); 

  mng_ptr_->num_files_ = prsm_ptrs.size() * 2 + PrsmUtil::getSpeciesIds(prsm_ptrs).size() +
      FastaUtil::countProteinNum(mng_ptr_->prsm_para_ptr_->getSearchDbFileName()) * 2;

  outputPrsms(prsm_ptrs);
  outputProteoforms(prsm_ptrs);
  outputProteins(prsm_ptrs);
  outputAllProteins(prsm_ptrs);
  outputAllPrsms(prsm_ptrs);
  outputFileList();
}

} /* namespace prot */
