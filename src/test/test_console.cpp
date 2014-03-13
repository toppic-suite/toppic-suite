#include <iostream>

#include "base/base_data.hpp"

#include "prsm/prsm_combine.hpp"
#include "prsm/prsm_selector.hpp"

#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_processor.hpp"

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"

#include "tdgf/evalue_processor.hpp"
#include "tdgf/tdgf_mng.hpp"

#include "console/output_selector.hpp"
#include "console/table_writer.hpp"

int main(int argc, char* argv[]) {
  try {
    std::cout << "Test msalign+ 0.1 " << std::endl;

    std::cout << "Zero ptm search " << std::endl;
    prot::ZeroPtmMngPtr zero_mng_ptr
        = prot::ZeroPtmMngPtr(new prot::ZeroPtmMng (std::string(argv[1])));
    zero_mng_ptr->search_db_file_name_ = argv[2];
    zero_mng_ptr->spectrum_file_name_ = argv[3];
    zero_mng_ptr->output_file_ext_ = "ZERO";
    prot::zeroPtmSearchProcess(zero_mng_ptr);

    std::cout << "Fast filter 0.1 " << std::endl;
    prot::PtmFastFilterMngPtr filter_mng_ptr = prot::PtmFastFilterMngPtr(new prot::PtmFastFilterMng(std::string(argv[1])));
		filter_mng_ptr->search_db_file_name_ = argv[2];
		filter_mng_ptr->spectrum_file_name_ = argv[3];
    filter_mng_ptr->output_file_ext_ = "FILTER";
    prot::PtmFastFilterProcessor filter_processor(filter_mng_ptr);
    filter_processor.process();


    std::cout << "Ptm alignment 0.1 " << std::endl;
		prot::PtmMngPtr ptm_mng_ptr = prot::PtmMngPtr(new prot::PtmMng(std::string(argv[1])));
		ptm_mng_ptr->search_db_file_name_ = argv[2];
		ptm_mng_ptr->spectrum_file_name_ = argv[3];
		ptm_mng_ptr->input_file_ext_ ="FILTER_COMBINED";
		ptm_mng_ptr->output_file_ext_="PTM";
    prot::PtmProcessor ptm_processor(ptm_mng_ptr);
		ptm_processor.process();

    std::cout << "Combine prsms " << std::endl;
		std::vector<std::string> input_exts ;
		input_exts.push_back("ZERO");
		input_exts.push_back("PTM");
    prot::PrSMCombine combine_processor(argv[2], argv[3], input_exts, "RAW_RESULT");
		combine_processor.process();

    std::cout << "EValueConsole 0.1 " << std::endl;
    prot::TdgfMngPtr tdgf_mng_ptr = prot::TdgfMngPtr(new prot::TdgfMng (std::string(argv[1])));
    tdgf_mng_ptr->search_db_file_name_ = argv[2];
    tdgf_mng_ptr->spectrum_file_name_ = argv[3];
    tdgf_mng_ptr->input_file_ext_ = "RAW_RESULT";
    tdgf_mng_ptr->output_file_ext_ = "EVALUE";
    prot::EValueProcessor processor(tdgf_mng_ptr);
    processor.init();
    /* compute E-value for a set of prsms each run */
    processor.process(false);


    std::cout << "Top selector 0.1 " << std::endl;
    prot::PrSMSelector selector(argv[2], argv[3], "EVALUE", "TOP", 1, "e");
		selector.process();

    std::cout << "Cutoff selector 0.1 " << std::endl;
    prot::OutputSelector output_selector(argv[2], argv[3], "TOP", "OUTPUT_RESULT", 
                                   "EVALUE", 0.01, 0.01, zero_mng_ptr->ppo_);
		output_selector.process();

    std::cout << "Table output 0.1" << std::endl;
    prot::TableWriter table_out(argv[2], argv[3], "OUTPUT_RESULT",
                          "OUTPUT_TABLE", zero_mng_ptr->ppo_);
		table_out.write();

  } catch (const char* e) {
    std::cout << "Exception " << e << std::endl;
  }
  return 0;
}
