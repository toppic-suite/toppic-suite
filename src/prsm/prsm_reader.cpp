#include <xercesc/framework/MemBufInputSource.hpp>

#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "prsm/prsm_reader.hpp"

namespace prot {

PrsmReader::PrsmReader(const std::string &file_name) {
  input_.open(file_name.c_str(), std::ios::in);
}

std::vector<std::string> PrsmReader::readOnePrsmLines() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    line = trim(line);
    if (line ==  "<prsm>") {
      line_list.push_back(line);
    }
    else if (line == "</prsm>") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    }
    else if (line == "") {
      continue;
    }
    else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

PrsmStrPtr PrsmReader::readOnePrsmStr() {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  if (prsm_str_vec.size() == 0) {
    return PrsmStrPtr(nullptr);
  }
  return PrsmStrPtr(new PrsmStr(prsm_str_vec));
}

PrsmPtr PrsmReader::readOnePrsm() {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  if (prsm_str_vec.size() == 0) {
    return PrsmPtr(nullptr);
  }
  std::string prsm_str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
  for (size_t i = 0; i < prsm_str_vec.size(); i++) {
    prsm_str += prsm_str_vec[i];
  }
  xercesc::MemBufInputSource prsm_buf(
      (const XMLByte*)prsm_str.c_str(), prsm_str.size(), "prsm_str (in memory)");

  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  PrsmPtr ptr;
  if(parser){
    XmlDOMDocument doc(parser, prsm_buf);
    xercesc::DOMElement* root = doc.getDocumentElement();
    //ptr = PrsmPtr(new Prsm(root));
  }
  //LOG_DEBUG("simple prsm spectrum id " << ptr->getSpectrumId() << " seq name " << ptr->getSeqName());
  return ptr;
}

/*
  xercesc::MemBufInputSource prsm_buf(
      (const XMLByte*)prsm_str.c_str(), prsm_str.size(), "prsm_str (in memory)");

  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  PrsmPtr ptr;
  if(parser){
    //XmlDOMDocument doc(parser, prsm_buf);
    //xercesc::DOMElement* root = doc.getDocumentElement();
    //ptr = PrsmPtr(new Prsm(root, proteo_ptrs));
  }
  return ptr;
}
*/

void PrsmReader::close() {
  input_.close();
}


} /* namespace prot */
