//
// Created by Abdul on 4/1/2024.
//

#ifndef TOPPIC_TOPDIA_SINGLE_PROCESS_HPP
#define TOPPIC_TOPDIA_SINGLE_PROCESS_HPP

#include <string>
#include <vector>
#include <map>

namespace toppic {

    namespace topdia_single_process {

        std::string geneArgumentStr(std::map<std::string, std::string> arguments,
                                    const std::string & prefix);

        int processOneFile(std::map<std::string, std::string> arguments,
                           const std::string &argument_str,
                           const std::string &spec_file_name, int frac_id);

        int process(std::map<std::string, std::string> arguments,
                    std::vector<std::string> spec_file_lst);

    }

}

#endif //TOPPIC_TOPDIA_SINGLE_PROCESS_HPP
