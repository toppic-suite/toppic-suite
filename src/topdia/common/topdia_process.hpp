//
// Created by Abdul on 4/1/2024.
//

#ifndef TOPPIC_TOPDIA_PROCESS_HPP
#define TOPPIC_TOPDIA_PROCESS_HPP

#include <string>
#include <vector>

#include "topdia/common/topdia_para.hpp"

namespace toppic {

namespace topdia_process {

int process(TopdiaParaPtr para_ptr,
            std::vector<std::string> spec_file_lst);

}

}

#endif //TOPPIC_TOPDIA_PROCESS_HPP
