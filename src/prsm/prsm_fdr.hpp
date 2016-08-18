#ifndef PROT_PRSM_PRSM_FDR_HPP_
#define PROT_PRSM_PRSM_FDR_HPP_

#include <map>
#include <string>

#include "base/string_util.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace prot {

class PrsmFdr {
 public:
  PrsmFdr(const std::string &db_file_name,
          const std::string &spec_file_name,
          const std::string &input_file_ext,
          const std::string &output_file_ext);
  void process();
 private:
  std::string db_file_name_;
  std::string spec_file_name_;
  std::string input_file_ext_;
  std::string output_file_ext_;

  void computeFdr(PrsmStrPtrVec &target, PrsmStrPtrVec &decoy);

  void computeProteoformFdr(PrsmStrPtr2D &target, PrsmStrPtr2D &decoy);
};
typedef std::shared_ptr<PrsmFdr> PrsmFdrPtr;

} /* namespace prot */

#endif /* PRSMFDR_HPP_ */
