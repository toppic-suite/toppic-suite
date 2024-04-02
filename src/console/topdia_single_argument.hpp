//
// Created by Abdul on 4/1/2024.
//

#ifndef TOPPIC_TOPDIA_SINGLE_ARGUMENT_HPP
#define TOPPIC_TOPDIA_SINGLE_ARGUMENT_HPP


#include <memory>
#include <vector>
#include <string>
#include <map>

#include <boost/program_options.hpp>

#include "common/xml/xml_dom_document.hpp"

namespace toppic {

    class Argument {
    public:
        Argument();

        bool parse(int argc, char* argv[]);

        std::map<std::string,std::string> getArguments(){ return arguments_;}

        std::vector<std::string> getSpecFileList() { return spec_file_list_;};
    private:
        void initArguments();

        void setArgumentsByConfigFile(const std::string &file_name);

        bool validateArguments();

        void showUsage(boost::program_options::options_description &desc);

    private:

        std::map<std::string,std::string> arguments_;

        std::vector<std::string> spec_file_list_;
    };

}  // namespace toppic

#endif //TOPPIC_TOPDIA_SINGLE_ARGUMENT_HPP
