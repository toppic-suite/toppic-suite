#include <xercesc/framework/MemBufInputSource.hpp>

#include "htslib/faidx.h"

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
    line = StringUtil::trim(line);
    //LOG_DEBUG("line " << line);
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

PrsmPtr PrsmReader::readOnePrsm(FastaIndexReaderPtr reader_ptr, 
                                const ModPtrVec fix_mod_list) {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  //LOG_DEBUG("prsm str vec size " << prsm_str_vec.size());
  if (prsm_str_vec.size() == 0) {
    return PrsmPtr(nullptr);
  }
  std::string prsm_str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
  for (size_t i = 0; i < prsm_str_vec.size(); i++) {
    prsm_str += prsm_str_vec[i];
  }
  //LOG_DEBUG("prsm str " << prsm_str);
  xercesc::MemBufInputSource prsm_buf(
      (const XMLByte*)prsm_str.c_str(), prsm_str.size(), "prsm_str (in memory)");

  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  PrsmPtr ptr;
  if(parser){
    XmlDOMDocument doc(parser, prsm_buf);
    xercesc::DOMElement* root = doc.getDocumentElement();
    ptr = PrsmPtr(new Prsm(root, reader_ptr, fix_mod_list));
  }
  //LOG_DEBUG("simple prsm spectrum id " << ptr->getSpectrumId() << " seq name " << ptr->getSeqName());
  return ptr;
}

void PrsmReader::close() {
  input_.close();
}

PrsmStrPtrVec readAllPrsmStrs(const std::string &input_file_name) {
  PrsmReader reader(input_file_name);
  PrsmStrPtrVec prsm_str_ptrs;
  PrsmStrPtr prsm_str_ptr = reader.readOnePrsmStr();
  while (prsm_str_ptr != nullptr) {
    prsm_str_ptrs.push_back(prsm_str_ptr);
    prsm_str_ptr = reader.readOnePrsmStr();
  }
  reader.close();
  return prsm_str_ptrs;
}

PrsmPtrVec readAllPrsms(const std::string &prsm_file_name, 
                        const std::string &db_file_name,
                        const ModPtrVec  &fix_mod_list) {
  FastaIndexReaderPtr fasta_reader_ptr(new FastaIndexReader(db_file_name));
  PrsmReader reader(prsm_file_name);
  PrsmPtrVec prsm_ptrs;
  PrsmPtr prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  while (prsm_ptr != nullptr) {
    prsm_ptrs.push_back(prsm_ptr);
    prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  }
  reader.close();
  return prsm_ptrs;
}


} /* namespace prot */
