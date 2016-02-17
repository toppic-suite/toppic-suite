#include "base/file_util.hpp"
#include "prsmview/transformer.hpp"
#include "prsmview/anno_view.hpp"

namespace prot
{

void translate(std::map<std::string,std::string> &arguments){

  std::string spectrum_file_name_ = arguments["spectrumFileName"];
  std::string xml_dir = FileUtil::basename(spectrum_file_name_) + "_xml";
  std::string html_dir = FileUtil::basename(spectrum_file_name_) + "_html";
  std::string exec_dir = arguments["executiveDir"];

  FileUtil::createFolder(html_dir + FileUtil::getFileSeparator() +"proteoforms");
  FileUtil::createFolder(html_dir + FileUtil::getFileSeparator() +"prsms");
  FileUtil::createFolder(html_dir + FileUtil::getFileSeparator() +"proteins");
  boost::filesystem::path from_path(exec_dir + FileUtil::getFileSeparator() + "toppic_resources" + FileUtil::getFileSeparator() + "web");
  boost::filesystem::path to_path(html_dir + FileUtil::getFileSeparator() + "resources");
  FileUtil::copyDir(from_path, to_path);

  //std::cout<<"trans start!XMLPlatformUtils::Initialize()"<<std::endl;
  xercesc::XMLPlatformUtils::Initialize();
  //std::cout<<"trans start! XalanTransformer::initialize()"<<std::endl;
  xalanc::XalanTransformer::initialize();
  //std::cout<<"trans start ! XalanTransformer"<<std::endl;
  xalanc::XalanTransformer theXanlanTransformer;

  std::string xml_file_list = xml_dir + FileUtil::getFileSeparator() + "files.xml";
  std::vector<std::vector<std::string>> anno_view = readViewXmlFiles(xml_file_list);

  for(size_t i = 0; i < anno_view.size(); i++) {
    std::cout << std::flush << "Processing " << i + 1 << " of " << anno_view.size() << " files.\r";
    const char* xml_in = anno_view[i][0].c_str();
    const char* xsl_in = anno_view[i][1].c_str();
    const char* xml_out = anno_view[i][2].c_str();

    LOG_DEBUG("xml in " << xml_in << " xsl in " << xsl_in << " xml out " << xml_out);

    theXanlanTransformer.transform(xml_in, xsl_in, xml_out);
  }

  std::cout << std::endl;
  xalanc::XalanTransformer::terminate();
  xercesc::XMLPlatformUtils::Terminate();
  xalanc::XalanTransformer::ICUCleanUp();
}

}
