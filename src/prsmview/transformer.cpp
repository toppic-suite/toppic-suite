#include "base/file_util.hpp"
#include "prsmview/transformer.hpp"
#include "prsmview/anno_view.hpp"

namespace prot
{

void translate(std::map<std::string,std::string> &arguments)
{
    std::string spectrum_file_name_ = arguments["spectrumFileName"];
    std::string xml_dir = basename(spectrum_file_name_) + "_xml";
    std::string html_dir = basename(spectrum_file_name_) + "_html";
    std::string exec_dir = arguments["executiveDir"];

    createFolder(html_dir + FILE_SEPARATOR +"proteoforms");
    createFolder(html_dir + FILE_SEPARATOR +"prsms");
    createFolder(html_dir + FILE_SEPARATOR +"proteins");
    boost::filesystem::path from_path(exec_dir + FILE_SEPARATOR + "resources");
    boost::filesystem::path to_path(html_dir + FILE_SEPARATOR + "resources");
    copyDir(from_path, to_path);



    //std::cout<<"trans start!XMLPlatformUtils::Initialize()"<<std::endl;
    xercesc::XMLPlatformUtils::Initialize();
    //std::cout<<"trans start! XalanTransformer::initialize()"<<std::endl;
    xalanc::XalanTransformer::initialize();
    //std::cout<<"trans start ! XalanTransformer"<<std::endl;
    xalanc::XalanTransformer theXanlanTransformer;

    std::string xml_file_list = xml_dir + FILE_SEPARATOR + "files.xml";
    std::vector<std::vector<std::string>> anno_view = readViewXmlFiles(xml_file_list);
    for(unsigned int i=0; i<anno_view.size(); i++)
    {
        //std::cout<<anno_view[i][0]<<std::endl;
        const char* xml_in = anno_view[i][0].c_str();
        const char* xsl_in = anno_view[i][1].c_str();
        const char* xml_out = anno_view[i][2].c_str();
        LOG_DEBUG("xml in " << xml_in << " xsl in " << xsl_in << " xml out " << xml_out);

        theXanlanTransformer.transform(xml_in,xsl_in,xml_out);
    }

    xalanc::XalanTransformer::terminate();
    xercesc::XMLPlatformUtils::Terminate();
    xalanc::XalanTransformer::ICUCleanUp();

}

}
