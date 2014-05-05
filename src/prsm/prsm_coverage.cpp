#include <iostream>
#include <fstream>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "prsm/prsm_coverage.hpp"

namespace prot {

PrsmCoverage::PrsmCoverage(PrsmParaPtr prsm_para_ptr,
                           std::string input_file_ext,
                           std::string output_file_ext) {
  prsm_para_ptr_ = prsm_para_ptr;
  input_file_ext_ = input_file_ext;
  output_file_ext_ = output_file_ext;
}

void PrsmCoverage::process(){

  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(prsm_para_ptr_->getSearchDbFileName(), 
                              prsm_para_ptr_->getFixModResiduePtrVec());

  LOG_DEBUG("protein data set loaded");
  std::string base_name = basename(prsm_para_ptr_->getSpectrumFileName());
  std::string input_file_name = base_name + "." + input_file_ext_;
  PrsmPtrVec prsms = readPrsm(input_file_name, raw_forms);
  LOG_DEBUG("read prsm complete ");
  addSpectrumPtrsToPrsms(prsms, prsm_para_ptr_);
  LOG_DEBUG("prsms loaded");

  std::string output_file_name = base_name+"."+output_file_ext_;
  std::ofstream file; 
  file.open(output_file_name.c_str());
  //write title
  file << "Data_file_name" << "\t"
      << "Prsm_ID" << "\t"
      << "Spectrum_ID"<< "\t"
      << "Activation_type" << "\t"
      << "Scan(s)" << "\t"
      << "#peaks"<< "\t"
      << "Charge" << "\t"
      << "Precursor_mass" << "\t"
      << "Adjusted_precursor_mass" << "\t"
      << "Protein_ID" << "\t"
      << "Species_ID" << "\t"
      << "Protein_name" << "\t"
      << "Protein_mass" << "\t"
      << "First_residue" << "\t"
      << "Last_residue" << "\t"
      << "Peptide" << "\t"
      << "#unexpected_modifications" << "\t"
      << "#matched_peaks" << "\t"
      << "#matched_fragment_ions" << "\t"
      << "P-Value" << "\t"
      << "E-Value" << "\t"
      << "One_Protein_probabilty"<< "\t"
      << "FDR" << "\t"

      << "N_term_ion_coverage" << "\t"
      << "C_term_ion_coverage" << "\t"
      << "Both_term_ion_coverage" << "\t"
      << "left_N_term_ion_coverage" << "\t"
      << "left_C_term_ion_coverage" << "\t"
      << "left_Both_term_ion_coverage" << "\t"
      << "middle_N_term_ion_coverage" << "\t"
      << "middle_C_term_ion_coverage" << "\t"
      << "middle_Both_term_ion_coverage" << "\t"
      << "right_N_term_ion_coverage" << "\t"
      << "right_C_term_ion_coverage" << "\t"
      << "right_Both_term_ion_coverage" << "\t"
      << std::endl;

  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  for(unsigned int i=0;i<prsms.size();i++){
    //std::cout << "proteom ptr " << prsms[i]->getProteoformPtr() << std::endl;
    //std::cout << "ms ptr " << prsms[i]->getRefineMs() << std::endl;

    PeakIonPairPtrVec pairs =  getPeakIonPairs (prsms[i]->getProteoformPtr(), 
                                                prsms[i]->getRefineMs(),
                                                min_mass);
    //std::cout << "get pair complete " << i  << std::endl;
    int len = prsms[i]->getProteoformPtr()->getResSeqPtr()->getLen() - 1;
    int begin = 1;
    int end = len - 1;
    double n_full_coverage = computePairConverage(pairs, begin, end, N_TERM_COVERAGE);
    double c_full_coverage = computePairConverage(pairs, begin, end, C_TERM_COVERAGE);
    double both_full_coverage = computePairConverage(pairs, begin, end, BOTH_TERM_COVERAGE);
    int one_third = len /3;
    end = one_third;
    double left_n_full_coverage = computePairConverage(pairs, begin, end, N_TERM_COVERAGE);
    double left_c_full_coverage = computePairConverage(pairs, begin, end, C_TERM_COVERAGE);
    double left_both_full_coverage = computePairConverage(pairs, begin, end, BOTH_TERM_COVERAGE);
    begin = one_third + 1;
    int two_thirds = len/3 * 2;
    end = two_thirds;
    double middle_n_full_coverage = computePairConverage(pairs, begin, end, N_TERM_COVERAGE);
    double middle_c_full_coverage = computePairConverage(pairs, begin, end, C_TERM_COVERAGE);
    double middle_both_full_coverage = computePairConverage(pairs, begin, end, BOTH_TERM_COVERAGE);
    begin = two_thirds + 1;
    end = len -1;
    double right_n_full_coverage = computePairConverage(pairs, begin, end, N_TERM_COVERAGE);
    double right_c_full_coverage = computePairConverage(pairs, begin, end, C_TERM_COVERAGE);
    double right_both_full_coverage = computePairConverage(pairs, begin, end, BOTH_TERM_COVERAGE);

    file << prsm_para_ptr_->getSpectrumFileName() << "\t"
        << prsms[i]->getId() << "\t"
        << prsms[i]->getSpectrumId()<< "\t"
        << prsms[i]->getDeconvMsPtr()->getHeaderPtr()->getActivationPtr()->getName()<< "\t"
        << prsms[i]->getSpectrumScan() << "\t"
        << prsms[i]->getDeconvMsPtr()->size()<< "\t"
        << prsms[i]->getDeconvMsPtr()->getHeaderPtr()->getPrecCharge() << "\t"
        << prsms[i]->getOriPrecMass()<< "\t"//"Precursor_mass"
        << prsms[i]->getAdjustedPrecMass() << "\t"
        << prsms[i]->getProteoformPtr()->getDbResSeqPtr()->getId() << "\t"
        << prsms[i]->getProteoformPtr()->getSpeciesId() << "\t"
        << prsms[i]->getProteoformPtr()->getDbResSeqPtr()->getName() << "\t"
        << prsms[i]->getProteoformPtr()->getDbResSeqPtr()->getSeqMass() << "\t"
        << prsms[i]->getProteoformPtr()->getStartPos() << "\t"
        << prsms[i]->getProteoformPtr()->getEndPos() << "\t"
        << prsms[i]->getProteoformPtr()->getProteinMatchSeq() << "\t"
        << prsms[i]->getProteoformPtr()->getUnexpectedChangeNum() << "\t"
        << prsms[i]->getMatchPeakNum() << "\t"
        << prsms[i]->getMatchFragNum() << "\t"
        << prsms[i]->getPValue() << "\t"
        << prsms[i]->getEValue() << "\t"
        << prsms[i]->getProbPtr()->getOneProtProb()<< "\t"
        << prsms[i]->getFdr() << "\t"

        << n_full_coverage << "\t"
        << c_full_coverage << "\t"
        << both_full_coverage << "\t"
        << left_n_full_coverage << "\t"
        << left_c_full_coverage << "\t"
        << left_both_full_coverage << "\t"
        << middle_n_full_coverage << "\t"
        << middle_c_full_coverage << "\t"
        << middle_both_full_coverage << "\t"
        << right_n_full_coverage << "\t"
        << right_c_full_coverage << "\t"
        << right_both_full_coverage << "\t"
        << std::endl;
    //std::cout << "print coverage complete " << std::endl;
  }
  file.close();
}

}
    
//    public static void outputCombineCoverage(List<Prsm> prsms, Properties prop, Map<Integer, Ms<DeconvPeak>> spectrums, AlignMng mng) throws Exception{
//        Collections.sort(prsms, new PrsmIdComparator());
//        int i = 0;
//        PrintWriter writer = new PrintWriter("msoutput/result_combine_coverage.txt");
//        PrsmTableWriter.combineHeading(writer);
//        int nSameForm = 0;
//        int nDiffFormCID = 0;
//        int nDiffFormETD = 0;
//        int nDiffProtCID = 0;
//        int nDiffProtETD = 0;
//        int nOneIdCID = 0;
//        int nOneIdETD = 0;
//        int count = 0;
//        while (i < prsms.size()-1) {
//            Prsm curPrsm = prsms.get(i);
//            Prsm nextPrsm = prsms.get(i+1);
//            int curScanNo = Integer.parseInt(curPrsm.getDeconvMs().getHeader().getScansString());
//            int nextScanNo = Integer.parseInt(nextPrsm.getDeconvMs().getHeader().getScansString());
//            EnumMsActivation nextType = nextPrsm.getDeconvMs().getHeader().getActivationType();
//
//            if (curScanNo + 1 == nextScanNo && nextType == EnumMsActivation.ETD) {
//                if (PrsmUtil.isCompatiablePtmSpecies(curPrsm, nextPrsm, mng.ppo )) {
//                    nSameForm++;
//                    PrsmTableWriter.writeCombineLine(writer, prop, curPrsm, nextPrsm);
//                }
//                else {
//                    if (curPrsm.getEValue() <= nextPrsm.getEValue()) {
//                        if (PrsmUtil.isSameProtein(curPrsm, nextPrsm)) {
//                            nDiffFormCID++;
//                            System.out.println("Same protein, different form: " 
//                                    + curPrsm.getDeconvMs().getHeader().getScansString() 
//                                    + " " + nextPrsm.getDeconvMs().getHeader().getScansString());
//                        }
//                        else {
//                            nDiffProtCID++;
//                        }
//                        Prsm newPrsm = PrsmReader.getOtherPrsm(prop, spectrums, curPrsm, mng);
//                        PrsmTableWriter.writeCombineLine(writer, prop, curPrsm, newPrsm);
//                    }
//                    else {
//                        if (PrsmUtil.isSameProtein(curPrsm, nextPrsm)) {
//                            nDiffFormETD++;
//                            System.out.println("Same protein, different form: " 
//                                    + curPrsm.getDeconvMs().getHeader().getScansString() 
//                                    + " " + nextPrsm.getDeconvMs().getHeader().getScansString());
//                        }
//                        else {
//                            nDiffProtETD++;
//                        }
//                        Prsm newPrsm = PrsmReader.getOtherPrsm(prop, spectrums, nextPrsm, mng);
//                        PrsmTableWriter.writeCombineLine(writer, prop, newPrsm, nextPrsm);
//                    }
//                }
//                count++;
//                i = i+2;
//            }
//            else {
//                
//                Prsm newPrsm = PrsmReader.getOtherPrsm(prop, spectrums, curPrsm, mng);
//                if (curPrsm.getDeconvMs().getHeader().getActivationType() == EnumMsActivation.ETD) {
//                    nOneIdETD++;
//                    PrsmTableWriter.writeCombineLine(writer, prop, newPrsm, curPrsm);
//                }
//                else {
//                    nOneIdCID++;
//                    PrsmTableWriter.writeCombineLine(writer, prop, curPrsm, newPrsm);
//                }
//                count++;
//                i++;
//            }
//            
//        }
//        System.out.println("Number of same protein form:" + nSameForm);
//        System.out.println("Number of same protein, but different protein form, CID's E-value < ETD's E-value:" + nDiffFormCID);
//        System.out.println("Number of same protein, but different protein form, CID's E-value > ETD's E-value:" + nDiffFormETD);
//        System.out.println("Number of different protein, CID's E-value < ETD's E-value:" + nDiffProtCID);
//        System.out.println("Number of different protein, CID's E-value > ETD's E-value:" + nDiffProtETD);
//        System.out.println("Number of one CID ID:" + nOneIdCID);
//        System.out.println("Number of one ETD ID:" + nOneIdETD);
//        writer.close();
//        System.out.println("Count " + count);
//    }

