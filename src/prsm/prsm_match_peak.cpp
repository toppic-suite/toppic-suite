

void XmlGenerator::outputPrsmMatchPeaks() {
  file_util::createFolder("match_peaks");
  std::string spectrum_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;

  PrsmReader prsm_reader(input_file_name);

  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_,
                                             mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getActivationPtr(),
                          mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getSkipList());
  SpectrumSetPtr spec_set_ptr;
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();

  size_t cnt = 0;
  while ((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0])!= nullptr) {
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec
            = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
        std::string file_name = mng_ptr_->xml_path_+ file_util::getFileSeparator() +
            "match_peaks" + file_util::getFileSeparator() + "scan_" +
            deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getScansStrting() + ".peak";

        geneAnnoPrsm(writer.getDoc(), prsm_ptr, mng_ptr_);

        cnt++;
        std::cout << std::flush << "Generating xml files - processing " << cnt << " PrSMs.\r";
        prsm_ptr = prsm_reader.readOnePrsm(fasta_reader_ptr_, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }
    }
  }
  std::cout << std::endl;

  prsm_reader.close();
  sp_reader.close();
}

