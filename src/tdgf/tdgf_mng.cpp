#include "tdgf/tdgf_mng.hpp"

namespace prot {

TdgfMng::TdgfMng(PrsmParaPtr prsm_para_ptr, 
                 int shift_num, double max_ptm_mass, bool use_gf,
                 bool variable_ptm,
                 const std::string &input_file_ext, 
                 const std::string &output_file_ext):
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext),
    prsm_para_ptr_(prsm_para_ptr),
    use_gf_(use_gf),
    variable_ptm_(variable_ptm),
    max_ptm_mass_(max_ptm_mass),
    unexpected_shift_num_(shift_num) {
    }

}
