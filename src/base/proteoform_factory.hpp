#ifndef PROT_BASE_PROTEOFORM_FACTORY_HPP_
#define PROT_BASE_PROTEOFORM_FACTORY_HPP_

#include "base/proteoform.hpp"

namespace prot {

class ProteoformFactory {
 public:
  /* get db proteoform */
  static ProteoformPtr geneDbProteoformPtr(FastaSeqPtr seq_ptr, ModPtrVec fix_mod_list);

  /* generate a proteoform with protein mod */
  static ProteoformPtr geneProtModProteoform(ProteoformPtr db_form_ptr,
                                             ProtModPtr prot_mod_ptr);

  static ProteoformPtrVec geneProtModProteoform(ProteoformPtr db_form_ptr,
                                                const ProtModPtrVec &prot_mod_ptrs);

  static ProteoformPtrVec2D gene2DProtModProteoform(const ProteoformPtrVec &db_form_ptrs,
                                                    const ProtModPtrVec &prot_mod_ptrs);
  /*
   * get subproteoform. local_start and local_end are relatively to
   * the start position in the original proteoform
   */
  static ProteoformPtr geneSubProteoform(ProteoformPtr proteoform_ptr,
                                         int local_start, int local_end);

  /* generate a proteoform vector with protein mod */
  static ProteoformPtrVec geneProtModProteoform(const ProteoformPtrVec &ori_forms,
                                                const ProtModPtrVec &prot_mods);

};

} /* namespace prot */

#endif /* PROTEOFORM_HPP_ */
