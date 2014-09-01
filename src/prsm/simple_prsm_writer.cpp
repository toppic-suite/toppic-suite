#include <iostream>
#include <algorithm>

#include "base/xml_dom_document.hpp"
#include "base/xml_dom.hpp"
#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/simple_prsm_writer.hpp"

namespace prot {

SimplePrsmWriter::SimplePrsmWriter(const std::string &file_name){
    file_.open(file_name.c_str());
    file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    file_ << "<simple_prsm_list>" << std::endl;
    XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
    doc_ = new XmlDOMDocument(impl->createDoc("simple_prsm_list"));
    serializer_ = impl->createSerializer();
}

SimplePrsmWriter::~SimplePrsmWriter(){
    serializer_->release();
    delete doc_;
}

void SimplePrsmWriter::close(){
  file_ << "</simple_prsm_list>" << std::endl;
  file_.close();
}

void SimplePrsmWriter::write(const SimplePrsmPtrVec &simple_prsm_ptrs){
  for(size_t i=0;i<simple_prsm_ptrs.size();i++){
    xercesc::DOMElement* element = simple_prsm_ptrs[i]->toXml(doc_);
    std::string str = writeToString(serializer_, element);
    writeToStreamByRemovingDoubleLF(file_, str);
    element->release();
  }
}

void combineBlock(const std::string &sp_file_name, int block_size, 
                  const std::string &output_file_ext, unsigned int result_num){
  SimplePrsmPtrVec2D match_ptrs;

  for(int i=0;i< block_size;i++){
    std::string block_file_name = basename(sp_file_name) + 
        "." + output_file_ext+"_"+std::to_string(i);
    match_ptrs.push_back(readSimplePrsms(block_file_name));
  }

  std::vector<int> pointers(block_size);

  MsAlignReader reader(sp_file_name);
  std::string output_file_name = basename(sp_file_name) 
      + "." + output_file_ext+"_COMBINED";
  SimplePrsmWriter writer(output_file_name);
  DeconvMsPtr deconv_ms_ptr;
  while((deconv_ms_ptr = reader.getNextMs()) != nullptr){
    SimplePrsmPtrVec selected_match_ptrs;
    for(size_t i = 0;i<match_ptrs.size();i++){
      for(size_t j = pointers[i]; j <match_ptrs[i].size(); j++){
        if(match_ptrs[i][j]->isSameSpectrum(deconv_ms_ptr->getHeaderPtr())){
          selected_match_ptrs.push_back(match_ptrs[i][j]);
        }
        if (match_ptrs[i][j]->isLargerSpectrumId(deconv_ms_ptr->getHeaderPtr())) {
          pointers[i] = j;
          break;
        }
      }
    }
    std::sort(selected_match_ptrs.begin(),selected_match_ptrs.end(),simplePrsmDown);
    SimplePrsmPtrVec result_match_ptrs;
    for(size_t i=0; i < selected_match_ptrs.size(); i++){
      if( i >= result_num) {
        break;
      }
      result_match_ptrs.push_back(selected_match_ptrs[i]);
    }

    writer.write(result_match_ptrs);
  }
  reader.close();
  writer.close();
}

} /* namespace prot */
