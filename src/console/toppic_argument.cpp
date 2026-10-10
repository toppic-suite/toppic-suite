// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "console/toppic_argument.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <set>
#include <sstream>
#include <string>
#include <utility>

#include "common/base/mod_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/str_util.hpp"
#include "common/util/version.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "filter/tag/tag_filter_mng.hpp"
#include "filter/tag/tag_index.hpp"
#include "para/sp_para.hpp"
#include "prsm/prsm_post_mass_match.hpp"
#include "seq/fasta_util.hpp"

namespace toppic {

namespace {

// A built-in PTM's mass as ptm_base.xml gives it: 15 significant digits
// write each mass of the file back exactly
std::string ptmMassStr(const PtmPtr& ptm_ptr) {
  std::ostringstream stream;
  stream << std::setprecision(15) << ptm_ptr->getMonoMass();
  return stream.str();
}

}  // namespace

ToppicArgument::ToppicArgument() { arguments_ = initArguments(); }

std::map<std::string, std::string> ToppicArgument::initArguments() {
  std::map<std::string, std::string> arguments;
  arguments["executiveDir"] = ".";
  arguments["resourceDir"] = "";
  arguments["oriDatabaseFileName"] = "";
  arguments["databaseFileName"] = "";
  arguments["spectrumFileName"] = "";
  arguments["databaseBlockSize"] = "60000000";
  arguments["maxFragmentLength"] = "500";
  arguments["minBlockNum"] = "20";
  arguments["activation"] = "FILE";
  arguments["fixedMod"] = "";
  arguments["allowProtMod"] = "NONE,NME,NME_ACETYLATION,M_ACETYLATION";
  arguments["allowProtType"] = "COMPLETE,PREFIX,SUFFIX,INTERNAL";
  arguments["nTermLabelMass"] = "0";
  arguments["useApproxSpectra"] = "false";
  arguments["variablePtmNum"] = "3";
  arguments["variablePtmFileName"] = "";
  // TopMG's identification hosted in toppic, behind the hidden --topmg
  // option: "disabled", "pipeline" (after toppic's own identification
  // stages) or "topmg-only" (the only database search: toppic's own
  // identification stages are skipped). Its options below take TopMG's own
  // defaults, so that the two agree when neither is set.
  arguments["topmgMode"] = "disabled";
  // the E-values of that stage: the generating function of toppic's own
  // stages, or with the hidden --topmg-use-mcmc TopMG's MCMC
  arguments["topmgUseMcmc"] = "false";
  arguments["topmgShiftNumber"] = "1";
  // the error tolerance (Da) of that stage's graph alignment
  // (GraphAlignMng::topmg_align_error_tolerance_), TopMG's own default
  arguments["topmgAlignErrorTolerance"] = "0.1";
  // the results the ASF filtering of that stage keeps per spectrum, TopMG's
  // own number; toppic's own stages keep filteringResultNumber
  arguments["asfFilteringResultNumber"] = "15";
  arguments["proteoGraphGap"] = "40";
  arguments["varPtmNumInGap"] = "5";
  arguments["useAsfDiag"] = "false";
  arguments["shiftNumber"] = "1";
  arguments["minShiftMass"] = "-500";
  arguments["maxShiftMass"] = "500";

  arguments["searchType"] = "TARGET";
  arguments["numOfTopPrsms"] = "1";
  arguments["filteringResultNumber"] = "20";
  arguments["precErrorTolerance"] = "10";
  // a precursor error tolerance in Dalton (hidden --prec-dalton-error) in
  // place of precErrorTolerance's ppm; empty: the ppm tolerance
  arguments["precDaltonError"] = "";
  arguments["fragErrorTolerance"] = "10";
  arguments["envCnnCutoff"] = "0.4";
  arguments["proteoformPpmError"] = "false";
  arguments["proteoformErrorTolerance"] = "1.2";
  arguments["cutoffSpectralType"] = "EVALUE";
  arguments["cutoffSpectralValue"] = "0.01";
  arguments["cutoffProteoformType"] = "EVALUE";
  arguments["cutoffProteoformValue"] = "0.01";
  arguments["cutoffProteinType"] = "EVALUE";
  arguments["cutoffProteinValue"] = "0.01";
  arguments["localPtmFileName"] = "";
  arguments["localThreshold"] = "0.15";
  // the localization of -B by the variable PTM search of stat/local_align
  // (the default since toppic_test experiment 140); hidden --local-original
  // selects the MIScore candidate search of stat/local (hidden --local-align
  // is accepted and selects the default)
  arguments["localAlign"] = "true";

  arguments["groupSpectrumNumber"] = "1";
  arguments["combinedOutputName"] = "";
  arguments["threadNumber"] = "1";
  arguments["useFeatureFile"] = "true";
  arguments["keepTempFiles"] = "false";
  arguments["keepDecoyResults"] = "false";
  arguments["useFdrGroups"] = "false";
  // the FDR ranks the PrSMs by a semi-supervised rescoring
  // (stat/rescore/prsm_rescore) instead of the E-value; on by default in a
  // target-decoy search (-d), off with the hidden --disable-rescore
  arguments["rescore"] = "true";
  // hidden --second-precursor: the spectra with two or more precursors are
  // searched again with the second one, without the first pass's matched
  // masses, and the two passes are reported together (prsm_second_prec)
  arguments["secondPrecursor"] = "false";
  // hidden --mass-calib: the two-pass m/z recalibration (ms/calib,
  // prsm/prsm_mass_calib): a zero-shift calibration search, the MS1 and MS/MS
  // error models fitted on its confident PrSMs, and the search on the
  // corrected spectra (<base>_calib_ms2.msalign) with the corrected copy of
  // TopFD's database (<base>_calib.sqlite) for the post mass matching and
  // the PrSM tables.
  arguments["massCalib"] = "false";
  arguments["combineResultOnly"] = "false";
  arguments["outputRawPrsmTable"] = "false";
  arguments["outputPrsmCoverage"] = "false";
  arguments["postMassMatch"] = "true";
  arguments["averageMass"] = "false";
  arguments["postMinPeakNum"] = "1";
  arguments["zeroShiftPrecOneDaError"] = "true";
  arguments["oneShiftPrecOneDaError"] = "false";
  arguments["varPtmPrecOneDaError"] = "false";
  arguments["fragMassExtension"] = "true";
  arguments["fragExtendMinMass"] = "5000";
  arguments["varPtmFilterCandNum"] = "";
  arguments["deamidSelector"] = "false";
  arguments["varPtmEvalueModResidues"] = "true";
  arguments["fragMassDedup"] = "false";
  arguments["multiSkipEvalue"] = "";
  arguments["multiTypeAlign"] = "true";
  // the second alignment of the shift searches at the refined precursor
  // (hidden --second-align): disabled, strict (the fragment tolerance on
  // every peak) or relax (the first alignment's tolerances)
  arguments["secondAlign"] = "true";
  // the zero-shift and variable PTM filters as one pass over the index
  // (hidden --separate-var-ptm-filter runs them one after the other;
  // toppic_test experiment 150)
  arguments["separateVarPtmFilter"] = "false";
  // the tolerance of the complementary (REVERSED) PRM peaks in the E-values'
  // generating function (hidden --evalue-complement-tolerance): relax (the
  // default), strict or separate (PrsmPara::useEvalueStrictComplement,
  // useEvalueSeparate; toppic_test experiment 146)
  arguments["evalueComplementTolerance"] = "relax";
  // the filters' mass list (hidden --filter-mass-list): merge, fragment or
  // single (SpPara::FilterMassListMode)
  arguments["filterMassList"] = "fragment";
  // the tag filter (hidden --tag-filter): the stages whose candidate lists
  // it joins (var, one, multi, topmg or all; empty = off), the tag length
  // and the proteins kept per spectrum (filter/tag; toppic_test experiment
  // 150)
  arguments["tagFilterStages"] = "";
  arguments["tagFilterLength"] = "4";
  arguments["tagFilterCandNum"] = "10";

  arguments["version"] = "";
  return arguments;
}

void ToppicArgument::outputArguments(
    std::ostream& output, const std::string& sep,
    std::map<std::string, std::string> arguments) {
  int gap = 46;
  output << "********************** Parameters **********************"
         << std::endl;
  output << std::setw(gap) << std::left << "Protein database file:" << sep
         << arguments["oriDatabaseFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Spectrum file:" << sep
         << arguments["spectrumFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Spectrum file masses:" << sep
         << (arguments["averageMass"] == "true" ? "Average" : "Monoisotopic")
         << std::endl;
  output << std::setw(gap) << std::left << "Number of combined spectra:" << sep
         << arguments["groupSpectrumNumber"] << std::endl;
  output << std::setw(gap) << std::left << "Fragmentation method:" << sep
         << arguments["activation"] << std::endl;
  output << std::setw(gap) << std::left << "Search type:" << sep
         << arguments["searchType"] << std::endl;

  if (arguments["fixedMod"] != "") {
    if (arguments["fixedMod"] == "C57") {
      output << std::setw(gap) << std::left << "Fixed modifications BEGIN"
             << std::endl;
      output << std::setw(gap) << std::left << "Carbamidomethylation" << sep
             << ptmMassStr(PtmBase::getPtmPtr_C57()) << sep << "C" << std::endl;
      output << std::setw(gap) << std::left << "Fixed modifications END"
             << std::endl;
    } else if (arguments["fixedMod"] == "C58") {
      output << std::setw(gap) << std::left << "Fixed modifications BEGIN"
             << std::endl;
      output << std::setw(gap) << std::left << "Carboxymethylation" << sep
             << ptmMassStr(PtmBase::getPtmPtr_C58()) << sep << "C" << std::endl;
      output << std::setw(gap) << std::left << "Fixed modifications END"
             << std::endl;
    } else {
      output << std::setw(gap) << std::left
             << "Fixed modifications file name:" << sep << arguments["fixedMod"]
             << std::endl;
      output << std::setw(gap) << std::left << "Fixed modifications BEGIN"
             << std::endl;
      std::vector<std::vector<std::string>> mod_data =
          mod_util::readModTxtForTsv(arguments["fixedMod"]);
      for (size_t i = 0; i < mod_data.size(); i++) {
        output << std::setw(gap) << std::left << mod_data[i][0] << sep
               << mod_data[i][1] << sep << mod_data[i][2] << std::endl;
      }
      output << std::setw(gap) << std::left << "Fixed modifications END"
             << std::endl;
    }
  }

  output << std::setw(gap) << std::left << "Allowed N-terminal forms:" << sep
         << arguments["allowProtMod"] << std::endl;

  output << std::setw(gap) << std::left << "Allowed proteoform types:" << sep
         << arguments["allowProtType"] << std::endl;

  // the FDR computation, with or without a variable PTM file
  if (arguments.at("useFdrGroups") == "true") {
    output << std::setw(gap) << std::left << "FDR computation:" << sep
           << "separate PrSM groups" << std::endl;
  }
  if (arguments.count("rescore") > 0 && arguments.at("rescore") == "true" &&
      arguments.at("searchType") == "TARGET+DECOY" &&
      arguments.at("useFdrGroups") != "true") {
    output << std::setw(gap) << std::left << "FDR computation:" << sep
           << "rescored PrSMs (linear SVM, 3-fold cross-validation)"
           << std::endl;
  }
  // printed only when set, so that the default parameter block is unchanged
  if (arguments.count("massCalib") > 0 &&
      arguments.at("massCalib") == "true") {
    output << std::setw(gap) << std::left << "Mass calibration:" << sep
           << "two-pass (zero-shift calibration search)" << std::endl;
  }

  if (arguments["variablePtmFileName"] != "" &&
      arguments["variablePtmNum"] != "0") {
    output << std::setw(gap) << std::left
           << "Maximum number of variable modifications:" << sep
           << arguments["variablePtmNum"] << std::endl;
    output << std::setw(gap) << std::left
           << "Use approximate spectra in protein filtering:" << sep
           << arguments["useApproxSpectra"] << std::endl;
    if (arguments.count("secondPrecursor") > 0 &&
        arguments.at("secondPrecursor") == "true") {
      output << std::setw(gap) << std::left << "Second precursor search:"
             << sep << "yes" << std::endl;
    }
    // printed only when set, so that the default parameter block is
    // unchanged
    if (arguments["varPtmFilterCandNum"] != "") {
      output << std::setw(gap) << std::left
             << "Variable PTM filter candidates per spectrum:" << sep
             << arguments["varPtmFilterCandNum"]
             << " (complete,prefix/suffix,internal)" << std::endl;
    }
    if (arguments["deamidSelector"] == "true") {
      output << std::setw(gap) << std::left
             << "Deamidation PrSM selection:" << sep << "true" << std::endl;
    }
    if (arguments["varPtmEvalueModResidues"] == "false") {
      output << std::setw(gap) << std::left
             << "Variable PTM E-values with modified residues:" << sep
             << "false" << std::endl;
    }
    output << std::setw(gap) << std::left
           << "Variable modifications file name:" << sep
           << arguments["variablePtmFileName"] << std::endl;
    output << std::setw(gap) << std::left << "Variable modifications BEGIN"
           << std::endl;
    std::vector<std::vector<std::string>> mod_data =
        mod_util::readModTxtForTsv(arguments["variablePtmFileName"]);
    for (size_t i = 0; i < mod_data.size(); i++) {
      output << std::setw(gap) << std::left << mod_data[i][0] << sep
             << mod_data[i][1] << sep << mod_data[i][2] << std::endl;
    }
    output << std::setw(gap) << std::left << "Variable modifications END"
           << std::endl;
  }

  output << std::setw(gap) << std::left
         << "Maximum number of unexpected modifications:" << sep
         << arguments["shiftNumber"] << std::endl;
  output << std::setw(gap) << std::left
         << "Maximum mass shift of modifications:" << sep
         << arguments["maxShiftMass"] << " Da" << std::endl;
  output << std::setw(gap) << std::left
         << "Minimum mass shift of modifications:" << sep
         << arguments["minShiftMass"] << " Da" << std::endl;
  output << std::setw(gap) << std::left << "Spectrum-level cutoff type:" << sep
         << arguments["cutoffSpectralType"] << std::endl;
  output << std::setw(gap) << std::left << "Spectrum-level cutoff value:" << sep
         << arguments["cutoffSpectralValue"] << std::endl;
  output << std::setw(gap) << std::left
         << "Proteoform-level cutoff type:" << sep
         << arguments["cutoffProteoformType"] << std::endl;
  output << std::setw(gap) << std::left
         << "Proteoform-level cutoff value:" << sep
         << arguments["cutoffProteoformValue"] << std::endl;
  output << std::setw(gap) << std::left << "Protein-level cutoff type:" << sep
         << arguments["cutoffProteinType"] << std::endl;
  output << std::setw(gap) << std::left << "Protein-level cutoff value:" << sep
         << arguments["cutoffProteinValue"] << std::endl;

  if (arguments["secondAlign"] == "false") {
    output << std::setw(gap) << std::left
           << "Second alignment at the refined precursor:" << sep << "false"
           << std::endl;
  }
  if (arguments["separateVarPtmFilter"] == "true") {
    output << std::setw(gap) << std::left
           << "Zero-shift and variable PTM filters:" << sep << "separate"
           << std::endl;
  }
  if (arguments["evalueComplementTolerance"] != "relax") {
    output << std::setw(gap) << std::left
           << "E-value tolerance of complementary ions:" << sep
           << arguments["evalueComplementTolerance"] << std::endl;
  }
  if (arguments["precDaltonError"] != "") {
    output << std::setw(gap) << std::left
           << "Error tolerance for precursor masses:" << sep
           << arguments["precDaltonError"] << " Da" << std::endl;
  } else {
    output << std::setw(gap) << std::left
           << "Error tolerance for precursor masses:" << sep
           << arguments["precErrorTolerance"] << " ppm" << std::endl;
  }
  output << std::setw(gap) << std::left
         << "Error tolerance for fragment masses:" << sep
         << arguments["fragErrorTolerance"] << " ppm" << std::endl;
  output << std::setw(gap) << std::left << "EnvCNN score cutoff:" << sep
         << arguments["envCnnCutoff"] << std::endl;
  output << std::setw(gap) << std::left
         << "Fragment mass deduplication:" << sep
         << (arguments["fragMassDedup"] == "true" ? "True" : "False")
         << std::endl;
  // printed only when it is not the default, so that the default parameter
  // block is unchanged
  if (arguments["filterMassList"] != "" &&
      arguments["filterMassList"] != "fragment") {
    output << std::setw(gap) << std::left << "Filter mass list:" << sep
           << arguments["filterMassList"] << std::endl;
  }
  // printed only when the tag filter is on
  if (arguments["tagFilterStages"] != "") {
    output << std::setw(gap) << std::left << "Tag filter stages:" << sep
           << arguments["tagFilterStages"] << std::endl;
    output << std::setw(gap) << std::left << "Tag filter tag length:" << sep
           << arguments["tagFilterLength"] << std::endl;
    output << std::setw(gap) << std::left
           << "Tag filter candidates per spectrum:" << sep
           << arguments["tagFilterCandNum"] << std::endl;
  }
  if (arguments["postMassMatch"] == "true") {
    output << std::setw(gap) << std::left << "Post mass matching:" << sep
           << "True" << std::endl;
    output << std::setw(gap) << std::left
           << "Post mass matching minimum isotopic peaks:" << sep
           << arguments["postMinPeakNum"] << std::endl;
  } else {
    output << std::setw(gap) << std::left << "Post mass matching:" << sep
           << "False" << std::endl;
  }
  if (arguments["proteoformPpmError"] == "false") {
    output << std::setw(gap) << std::left
           << "Error tolerance for identifying PrSM clusters:" << sep
           << arguments["proteoformErrorTolerance"] << " Da" << std::endl;
  } else {
    output << std::setw(gap) << std::left
           << "Error tolerance for identifying PrSM clusters:" << sep
           << arguments["proteoformErrorTolerance"] << " PPM" << std::endl;
  }

  if (arguments["useFeatureFile"] == "true") {
    output << std::setw(gap) << std::left << "Use TopFD feature file:" << sep
           << "True" << std::endl;
  } else {
    output << std::setw(gap) << std::left << "Use TopFD feature file:" << sep
           << "False" << std::endl;
  }

  output << std::setw(gap) << std::left << "E-value computation:" << sep
         << "Generating function" << std::endl;

  if (arguments["localPtmFileName"] != "") {
    output << std::setw(gap) << std::left << "Localization with MIScore:" << sep
           << "True" << std::endl;
    output << std::setw(gap) << std::left << "MIScore threshold:" << sep
           << arguments["localThreshold"] << std::endl;
    output << std::setw(gap) << std::left << "Localization method:" << sep
           << (arguments["localAlign"] == "true"
                   ? "Variable PTM search (local_align)"
                   : "MIScore candidate search (local, --local-original)")
           << std::endl;
    output << std::setw(gap) << std::left
           << "Localization modification file name:" << sep
           << arguments["localPtmFileName"] << std::endl;
    output << std::setw(gap) << std::left << "Localization modifications BEGIN"
           << std::endl;
    std::vector<std::vector<std::string>> local_mod_data =
        mod_util::readModTxtForTsv(arguments["localPtmFileName"]);
    for (size_t i = 0; i < local_mod_data.size(); i++) {
      output << std::setw(gap) << std::left << local_mod_data[i][0] << sep
             << local_mod_data[i][1] << sep << local_mod_data[i][2]
             << std::endl;
    }
    output << std::setw(gap) << std::left << "Localization modifications END"
           << std::endl;
  } else {
    output << std::setw(gap) << std::left << "Localization with MIScore:" << sep
           << "False" << std::endl;
  }

  output << std::setw(gap) << std::left << "Thread number: " << sep
         << arguments["threadNumber"] << std::endl;

  output << std::setw(gap) << std::left << "Executable file directory:" << sep
         << arguments["executiveDir"] << std::endl;
  output << std::setw(gap) << std::left << "Start time:" << sep
         << arguments["startTime"] << std::endl;
  if (arguments["endTime"] != "") {
    output << std::setw(gap) << std::left << "End time:" << sep
           << arguments["endTime"] << std::endl;
  }
  output << std::setw(gap) << std::left << "Version:" << sep
         << arguments["version"] << std::endl;
  output << "********************** Parameters **********************"
         << std::endl;
}

std::string ToppicArgument::outputTsvArguments(
    std::map<std::string, std::string> arguments) {
  std::stringstream output;
  outputArguments(output, "\t", arguments);
  return output.str();
}

void ToppicArgument::showUsage(
    boost::program_options::options_description& desc) {
  std::cout << "Usage: toppic [options] database-file-name spectrum-file-name"
            << std::endl;
  std::cout << desc << std::endl;
  std::cout << "Version: " << Version::getVersion() << std::endl;
}

bool ToppicArgument::parse(int argc, char* argv[]) {
  std::string database_file_name = "";
  std::string argument_file_name = "";
  std::string activation = "";
  std::string fixed_mod = "";
  std::string allow_mod = "";
  std::string allow_prot_type = "";
  std::string n_term_label_mass = "";
  std::string variable_ptm_num = "";
  std::string variable_ptm_file_name = "";
  std::string shift_num = "";
  std::string topmg_mode = "";
  std::string topmg_shift_num = "";
  std::string topmg_align_error_tole = "";
  std::string asf_filtering_result_num = "";
  std::string min_shift_mass = "";
  std::string max_shift_mass = "";
  std::string frag_error_tole = "";
  std::string prec_error_tole = "";
  std::string env_cnn_cutoff = "";
  std::string form_error_tole = "";
  std::string cutoff_spectral_type = "";
  std::string cutoff_spectral_value = "";
  std::string cutoff_proteoform_type = "";
  std::string cutoff_proteoform_value = "";
  std::string cutoff_protein_type = "";
  std::string cutoff_protein_value = "";
  std::string local_ptm_file_name = "";
  std::string local_threshold = "";
  std::string group_num = "";
  std::string combined_output_name = "";
  std::string filtering_result_num = "";
  std::string thread_number = "";
  std::string top_prsm_number = "";
  std::string post_min_peak_num = "";
  std::string frag_extend_min_mass = "";
  std::string prec_dalton_error = "";
  std::string multi_skip_evalue = "";
  std::string var_ptm_filter_cand_num = "";
  std::string filter_mass_list = "";
  std::string evalue_complement_tolerance = "";
  std::string tag_filter_stages = "";
  std::string tag_filter_length = "";
  std::string tag_filter_cand_num = "";

  /** Define and parse the program options*/
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options()("help,h", "Print the help message.")(
        "activation,a", po::value<std::string>(&activation),
        "<CID|ETD|ETCID|ETHCD|HCD|UVPD|FILE>. Fragmentation method of MS/MS "
        "spectra. When FILE is used, fragmentation methods of spectra are "
        "given in the input spectral data file. ETCID (ETciD) and ETHCD "
        "(EThcD) are searched with the ion types of ETD (c and z-dot). Default "
        "value: FILE.")(
        "fixed-mod,f", po::value<std::string>(&fixed_mod),
        "<C57|C58|a fixed modification file>. Fixed modifications. Three "
        "available options: C57, C58, or the name of a text file containing "
        "the information of fixed modifications. When C57 is selected, "
        "carbamidomethylation on cysteine is the only fixed modification. When "
        "C58 is selected, carboxymethylation on cysteine is the only fixed "
        "modification.")(
        "n-terminal-form,n", po::value<std::string>(&allow_mod),
        "<a list of allowed N-terminal forms>. N-terminal forms of proteins. "
        "Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, "
        "and M_ACETYLATION. NONE stands for no modifications, NME for "
        "N-terminal methionine excision, NME_ACETYLATION for N-terminal "
        "acetylation after the initiator methionine is removed, and "
        "M_ACETYLATION for N-terminal methionine acetylation. When multiple "
        "forms are allowed, they are separated by commas. Default value: "
        "NONE,NME,NME_ACETYLATION,M_ACETYLATION.")(
        "n-terminal-label,N", po::value<std::string>(&n_term_label_mass),
        "<a non-negative number>. Mass of a label on the N-terminus of every "
        "proteoform, for example 229.162932 for TMT. The label is added to "
        "the N-terminal fragment ions and removed from the precursor mass; "
        "the reported proteoforms and adjusted precursor masses are "
        "label-free, and the result tables give the label in the column "
        "N-terminal label mass. Default value: 0 Dalton.")(
        "proteoform-type,R", po::value<std::string>(&allow_prot_type),
        "<a list of allowed types of proteoforms>. Four types of proteoforms "
        "can be selected: COMPLETE, PREFIX, SUFFIX, and INTERNAL. COMPLETE "
        "stands for proteoforms that cover the complete sequence of a protein "
        "(NME is allowed), PREFIX for proteoforms that cover the N-terminal "
        "part of a protein, SUFFIX for proteoforms that cover the C-terminal "
        "part of a protein, and INTERNAL for proteoforms that cover an "
        "internal part of a protein. When multiple types are allowed, they are "
        "separated by commas. Default value: COMPLETE,PREFIX,SUFFIX,INTERNAL.")(
        "num-shift,s", po::value<std::string>(&shift_num),
        "<0|1|2>. Maximum number of unexpected modifications in a "
        "proteoform-spectrum-match. Default value: 1.")(
        "min-shift,m", po::value<std::string>(&min_shift_mass),
        "Minimum value of the mass shift of an unexpected modification. "
        "Default value: -500 Dalton.")(
        "max-shift,M", po::value<std::string>(&max_shift_mass),
        "Maximum value of the mass shift of an unexpected modification. "
        "Default value: 500 Dalton.")
        ("variable-ptm-num,S", po::value<std::string>(&variable_ptm_num),
         "Maximum number of variable modifications in a "
         "proteoform-spectrum-match. Default value: 3.")(
            "variable-ptm-file-name,b",
            po::value<std::string>(&variable_ptm_file_name),
            "<a variable modification file>. Specify a text file containing "
            "the information of variable modifications.")(
            "decoy,d",
            "Use a shuffled decoy protein database to estimate false discovery "
            "rates.")("frag-error-tolerance,e",
                      po::value<std::string>(&frag_error_tole),
                      "<a positive integer>. Error tolerance for fragment "
                      "masses in PPM. Default value: 10.")(
            "prec-error-tolerance,E", po::value<std::string>(&prec_error_tole),
            "<a positive integer>. Error tolerance for precursor masses in "
            "PPM. Default value: 10.")(
            "filter-by-env-cnn,F", po::value<std::string>(&env_cnn_cutoff),
            "<a number between 0 and 1>. Remove the deconvoluted masses whose "
            "EnvCNN score is below the cutoff before spectral filtering, "
            "database search and E-value computation. Default value: 0.4.")(
            "proteoform-error-tolerance,p",
            po::value<std::string>(&form_error_tole),
            "<a positive number>. Error tolerance for identifying PrSM "
            "clusters. Default value: 1.2 Dalton.")(
            "spectrum-cutoff-type,t",
            po::value<std::string>(&cutoff_spectral_type),
            "<EVALUE|FDR>. Spectrum-level cutoff type for filtering identified "
            "proteoform-spectrum-matches. Default value: EVALUE.")(
            "spectrum-cutoff-value,v",
            po::value<std::string>(&cutoff_spectral_value),
            "<a positive number>. Spectrum-level cutoff value for filtering "
            "identified proteoform-spectrum-matches. Default value: 0.01.")(
            "proteoform-cutoff-type,T",
            po::value<std::string>(&cutoff_proteoform_type),
            "<EVALUE|FDR>. Proteoform-level cutoff type for filtering "
            "identified proteoform-spectrum-matches. Default value: EVALUE.")(
            "proteoform-cutoff-value,V",
            po::value<std::string>(&cutoff_proteoform_value),
            "<a positive number>. Proteoform-level cutoff value for filtering "
            "identified proteoform-spectrum-matches. Default value: 0.01.")(
            "protein-cutoff-type,y",
            po::value<std::string>(&cutoff_protein_type),
            "<EVALUE|FDR>. Protein-level cutoff type for filtering identified "
            "proteoform-spectrum-matches. With FDR, a protein is represented "
            "by its best proteoform (the lowest E-value). Default value: "
            "EVALUE.")(
            "protein-cutoff-value,Y",
            po::value<std::string>(&cutoff_protein_value),
            "<a positive number>. Protein-level cutoff value for filtering "
            "identified proteoform-spectrum-matches. Default value: 0.01.")(
            "approximate-spectra,A",
            "Use approximate spectra to increase the sensitivity in protein "
            "filtering. Default value: false.")(
            "local-ptm-file-name,B",
            po::value<std::string>(&local_ptm_file_name),
            "<a common modification file>. Specify a text file containing the "
            "information of common modifications for the characterization of "
            "unexpected shifts.")(
            "miscore-threshold,H", po::value<std::string>(&local_threshold),
            "<a positive number between 0 and 1>. Score threshold "
            "(modification identification score) for filtering results of "
            "modification characterization. Default value: 0.15.")(
            "thread-number,u", po::value<std::string>(&thread_number),
            "<a positive integer>. Number of threads used in the computation. "
            "Default value: 1.")(
            "num-combined-spectra,r", po::value<std::string>(&group_num),
            "<a positive integer>. Number of combined spectra. The parameter "
            "is set to 2 (or 3) for combining spectral pairs (or triplets) "
            "generated by the alternating fragmentation mode. Default value: "
            "1")("combined-file-name,c",
                 po::value<std::string>(&combined_output_name),
                 "Specify a file name for the combined spectrum data file and "
                 "analysis results.")(
            "no-topfd-feature,x",
            "No TopFD feature file for proteoform identification.")(
            "keep-temp-files,k", "Keep intermediate files.")(
            "keep-decoy-ids,K", "Keep decoy identifications.")(
            "disable-post-match,D",
            "Disable post mass matching. By default, after the search, the "
            "theoretical fragment masses of each PrSM that no deconvoluted "
            "mass matched are matched against the centroided MS/MS peaks in "
            "the TopFD SQLite database (run topfd without -N) with the method "
            "of MSPathFinderT, and the matched masses are added to the "
            "spectra (written to <name>_post_ms2.msalign).")(
            "post-min-peak-num,I", po::value<std::string>(&post_min_peak_num),
            "<a positive integer>. Minimum number of observed isotopic peaks "
            "of a fragment mass matched by the post mass matching. Default "
            "value: 1.")(
            "average-mass,g",
            "The spectrum files report average masses of precursor and "
            "fragment ions (topfd --average-mass). Before the search, each "
            "file is converted to <name>_mono_ms2.msalign with monoisotopic "
            "masses, which is then searched.");

    po::options_description desc("Options");

    desc.add_options()("help,h", "")("activation,a",
                                     po::value<std::string>(&activation), "")(
        "fixed-mod,f", po::value<std::string>(&fixed_mod), "")(
        "n-terminal-form,n", po::value<std::string>(&allow_mod), "")(
        "n-terminal-label,N", po::value<std::string>(&n_term_label_mass), "")(
        "num-shift,s", po::value<std::string>(&shift_num), "")(
        "min-shift,m", po::value<std::string>(&min_shift_mass), "")(
        "max-shift,M", po::value<std::string>(&max_shift_mass), "")(
        "variable-ptm-num,S", po::value<std::string>(&variable_ptm_num), "")(
        "variable-ptm-file-name,b",
        po::value<std::string>(&variable_ptm_file_name),
        "")("approximate-spectra,A", "")("decoy,d", "")(
        "frag-error-tolerance,e", po::value<std::string>(&frag_error_tole), "")(
        "prec-error-tolerance,E", po::value<std::string>(&prec_error_tole), "")(
        "prec-dalton-error", po::value<std::string>(&prec_dalton_error), "")(
        "filter-by-env-cnn,F", po::value<std::string>(&env_cnn_cutoff), "")(
        "proteoform-ppm-error,P", "")("proteoform-error-tolerance,p",
                                      po::value<std::string>(&form_error_tole),
                                      "")(
        "proteoform-type,R", po::value<std::string>(&allow_prot_type), "")(
        "spectrum-cutoff-type,t", po::value<std::string>(&cutoff_spectral_type),
        "")("spectrum-cutoff-value,v",
            po::value<std::string>(&cutoff_spectral_value),
            "")("proteoform-cutoff-type,T",
                po::value<std::string>(&cutoff_proteoform_type),
                "")("proteoform-cutoff-value,V",
                    po::value<std::string>(&cutoff_proteoform_value),
                    "")("protein-cutoff-type,y",
                        po::value<std::string>(&cutoff_protein_type), "")(
        "protein-cutoff-value,Y", po::value<std::string>(&cutoff_protein_value),
        "")("local-ptm-file-name,B",
            po::value<std::string>(&local_ptm_file_name), "")(
        "miscore-threshold,H", po::value<std::string>(&local_threshold), "")(
        "num-combined-spectra,r", po::value<std::string>(&group_num), "")(
        "combined-file-name,c", po::value<std::string>(&combined_output_name),
        "")("thread-number,u", po::value<std::string>(&thread_number), "")(
        "top-prsm-number,U", po::value<std::string>(&top_prsm_number), "")(
        "no-topfd-feature,x", "")("keep-temp-files,k", "")(
        "keep-decoy-ids,K", "")("combine-result-only,C", "")(
        "output-raw-prsm-table,o", "")("output-prsm-coverage,O", "")(
        "disable-post-match,D", "")(
        "post-min-peak-num,I", po::value<std::string>(&post_min_peak_num), "")(
        "average-mass,g", "")("disable-zero-shift-prec-error", "")(
        "enable-one-shift-prec-error", "")("enable-var-ptm-prec-error", "")(
        "fdr-groups", "")("disable-rescore", "")("second-precursor", "")(
        "mass-calib", "")(
        "local-align", "")(
        "local-original", "")(
        "disable-frag-mass-extension", "")(
        "frag-mass-dedup", "")(
        "frag-extend-min-mass", po::value<std::string>(&frag_extend_min_mass),
        "")("var-ptm-filter-cand-num",
            po::value<std::string>(&var_ptm_filter_cand_num), "")(
        "deamid-selector", "")("disable-var-ptm-evalue-mod-residues", "")(
        "multi-skip-evalue", po::value<std::string>(&multi_skip_evalue),
            "")("disable-multi-type-align", "")(
        "disable-second-align", "")("separate-var-ptm-filter", "")(
        "evalue-complement-tolerance",
        po::value<std::string>(&evalue_complement_tolerance), "")(
        "filter-mass-list", po::value<std::string>(&filter_mass_list), "")(
        "tag-filter", po::value<std::string>(&tag_filter_stages), "")(
        "tag-filter-length", po::value<std::string>(&tag_filter_length), "")(
        "tag-filter-cand-num", po::value<std::string>(&tag_filter_cand_num),
        "")("topmg", po::value<std::string>(&topmg_mode), "")(
        "topmg-use-mcmc", "")(
        "topmg-num-shift", po::value<std::string>(&topmg_shift_num), "")(
        "topmg-align-error-tolerance",
        po::value<std::string>(&topmg_align_error_tole), "")(
        "asf-filtering-result-number",
        po::value<std::string>(&asf_filtering_result_num), "")(
        "filtering-result-number",
        po::value<std::string>(&filtering_result_num),
        "Filtering result number. Default value: 20.")(
        "database-file-name",
        po::value<std::string>(&database_file_name)->required(),
        "Database file name with its path.")(
        "spectrum-file-name",
        po::value<std::vector<std::string>>()->multitoken()->required(),
        "Spectrum file name with its path.");

    po::positional_options_description positional_options;
    positional_options.add("database-file-name", 1);
    positional_options.add("spectrum-file-name", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv)
                    .options(desc)
                    .positional(positional_options)
                    .run(),
                vm);
      if (vm.count("help")) {
        showUsage(display_desc);
        return false;
      }
      po::notify(vm);
      // throws on error, so do after help in case there are any problems
    } catch (boost::program_options::required_option& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    } catch (boost::program_options::error& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }
    std::string argv_0(argv[0]);
    if (vm.count("full-binary-path")) {
      arguments_["executiveDir"] = argv[0];
    } else {
      arguments_["executiveDir"] = file_util::getExecutiveDir(argv_0);
      if (file_util::checkSpace(arguments_["executiveDir"])) {
        LOG_ERROR("Current directory "
                  << arguments_["executiveDir"]
                  << " contains space and will cause errors in the program!");
        exit(EXIT_FAILURE);
      }
    }
    LOG_DEBUG("Executive Dir " << arguments_["executiveDir"]);

    arguments_["resourceDir"] =
        file_util::getResourceDir(arguments_["executiveDir"]);

    arguments_["oriDatabaseFileName"] = database_file_name;

    if (vm.count("spectrum-file-name")) {
      spec_file_list_ = vm["spectrum-file-name"].as<std::vector<std::string>>();
    }

    if (vm.count("activation")) {
      arguments_["activation"] = activation;
    }

    if (vm.count("fixed-mod")) {
      arguments_["fixedMod"] = fixed_mod;
    }

    if (vm.count("n-terminal-form")) {
      arguments_["allowProtMod"] = allow_mod;
    }

    if (vm.count("proteoform-type")) {
      arguments_["allowProtType"] = allow_prot_type;
    }

    if (vm.count("n-terminal-label")) {
      arguments_["nTermLabelMass"] = n_term_label_mass;
    }

    if (vm.count("variable-ptm-num")) {
      arguments_["variablePtmNum"] = variable_ptm_num;
    }

    if (vm.count("variable-ptm-file-name")) {
      arguments_["variablePtmFileName"] = variable_ptm_file_name;
    }

    if (vm.count("approximate-spectra")) {
      arguments_["useApproxSpectra"] = "true";
    }

    if (vm.count("num-shift")) {
      arguments_["shiftNumber"] = shift_num;
    }

    if (vm.count("min-shift")) {
      arguments_["minShiftMass"] = min_shift_mass;
    }

    if (vm.count("max-shift")) {
      arguments_["maxShiftMass"] = max_shift_mass;
    }

    if (vm.count("decoy")) {
      arguments_["searchType"] = "TARGET+DECOY";
    }

    if (vm.count("frag-error-tolerance")) {
      arguments_["fragErrorTolerance"] = frag_error_tole;
    }

    if (vm.count("prec-error-tolerance")) {
      arguments_["precErrorTolerance"] = prec_error_tole;
    }

    // hidden: a precursor error tolerance in Dalton, used for every
    // precursor tolerance in place of -E's ppm (issue #563, toppic_test
    // experiment 133)
    if (vm.count("prec-dalton-error")) {
      arguments_["precDaltonError"] = prec_dalton_error;
    }

    if (vm.count("filter-by-env-cnn")) {
      arguments_["envCnnCutoff"] = env_cnn_cutoff;
    }

    if (vm.count("proteoform-ppm-error")) {
      arguments_["proteoformPpmError"] = "true";
    }

    if (vm.count("proteoform-error-tolerance")) {
      arguments_["proteoformErrorTolerance"] = form_error_tole;
    }

    if (vm.count("spectrum-cutoff-type")) {
      arguments_["cutoffSpectralType"] = cutoff_spectral_type;
    }

    if (vm.count("spectrum-cutoff-value")) {
      arguments_["cutoffSpectralValue"] = cutoff_spectral_value;
    }

    if (vm.count("proteoform-cutoff-type")) {
      arguments_["cutoffProteoformType"] = cutoff_proteoform_type;
    }

    if (vm.count("proteoform-cutoff-value")) {
      arguments_["cutoffProteoformValue"] = cutoff_proteoform_value;
    }

    if (vm.count("protein-cutoff-type")) {
      arguments_["cutoffProteinType"] = cutoff_protein_type;
    }

    if (vm.count("protein-cutoff-value")) {
      arguments_["cutoffProteinValue"] = cutoff_protein_value;
    }

    if (vm.count("local-ptm-file-name")) {
      arguments_["localPtmFileName"] = local_ptm_file_name;
    }

    if (vm.count("miscore-threshold")) {
      arguments_["localThreshold"] = local_threshold;
    }

    if (vm.count("num-combined-spectra")) {
      arguments_["groupSpectrumNumber"] = group_num;
    }

    if (vm.count("combined-file-name")) {
      arguments_["combinedOutputName"] = combined_output_name;
    }

    if (vm.count("thread-number")) {
      arguments_["threadNumber"] = thread_number;
    }

    if (vm.count("top-prsm-number")) {
      arguments_["numOfTopPrsms"] = top_prsm_number;
    }

    if (vm.count("no-topfd-feature")) {
      arguments_["useFeatureFile"] = "false";
    }

    if (vm.count("keep-temp-files")) {
      arguments_["keepTempFiles"] = "true";
    }
    if (vm.count("keep-decoy-ids")) {
      arguments_["keepDecoyResults"] = "true";
    }
    if (vm.count("fdr-groups")) {
      arguments_["useFdrGroups"] = "true";
    }
    if (vm.count("disable-rescore")) {
      arguments_["rescore"] = "false";
    }
    if (vm.count("second-precursor")) {
      arguments_["secondPrecursor"] = "true";
    }
    if (vm.count("mass-calib")) {
      arguments_["massCalib"] = "true";
    }
    // hidden: by default, the multiple-unexpected-shift search computes
    // COMPLETE/PREFIX/SUFFIX/INTERNAL with one fused dynamic program instead
    // of one PsAlign pass per type (PrsmPara::useMultiTypeAlign);
    // --disable-multi-type-align falls back to the original one-pass-per-type
    // PsAlign search.
    if (vm.count("disable-multi-type-align")) {
      arguments_["multiTypeAlign"] = "false";
    }

    // hidden: the one-shift and multiple-shift searches repeat their
    // alignment at the refined precursor of the first alignment's result,
    // with the fragment tolerance on every peak, and keep that precursor
    // (search/diag/second_align_util, issue #614, toppic_test experiment
    // 144); --disable-second-align keeps the first alignment's PrSM.
    if (vm.count("disable-second-align")) {
      arguments_["secondAlign"] = "false";
    }

    // hidden: the zero-shift and variable PTM filters one after the other
    // instead of one pass over the index files (toppic_test experiment 150)
    if (vm.count("separate-var-ptm-filter")) {
      arguments_["separateVarPtmFilter"] = "true";
    }

    // hidden: the tolerance of the complementary (REVERSED, P - m) PRM
    // peaks in the generating function of the E-values: relax (the
    // default, the precursor's error added), strict (the fragment
    // tolerance on the fragment mass) or separate (one generating function
    // per PrSM on PRM spectra at its own precursor, every peak strict;
    // issue #616, toppic_test experiment 146, report 0146: strict and
    // separate move the E-value scale without more identifications)
    if (vm.count("evalue-complement-tolerance")) {
      arguments_["evalueComplementTolerance"] = evalue_complement_tolerance;
    }

    // hidden: TopMG's identification - the ASF filtering with the variable
    // modification file, the graph alignment, its post-processing, the
    // E-value and the top PrSM selection. "pipeline" runs it after toppic's
    // own identification stages, on the spectra they left unidentified;
    // "topmg-only" runs it alone, on every spectrum, and skips toppic's own
    // identification stages (zero, variable-PTM, one and multiple unexpected
    // shift searches); "disabled", the default, does not run it.
    if (vm.count("topmg")) {
      arguments_["topmgMode"] = topmg_mode;
    }
    // hidden: compute the E-values of that stage with TopMG's MCMC instead
    // of the generating function (TDGF)
    if (vm.count("topmg-use-mcmc")) {
      arguments_["topmgUseMcmc"] = "true";
    }
    // hidden: the unexpected shifts of that graph alignment, TopMG's own
    // number, which --num-shift does not change
    if (vm.count("topmg-num-shift")) {
      arguments_["topmgShiftNumber"] = topmg_shift_num;
    }
    // hidden: the error tolerance (Da) of that graph alignment, TopMG's
    // hidden --align-error-tolerance
    if (vm.count("topmg-align-error-tolerance")) {
      arguments_["topmgAlignErrorTolerance"] = topmg_align_error_tole;
    }
    // hidden: the ASF filtering result number of that stage, TopMG's own,
    // which --filtering-result-number does not change
    if (vm.count("asf-filtering-result-number")) {
      arguments_["asfFilteringResultNumber"] = asf_filtering_result_num;
    }

    if (vm.count("combine-result-only")) {
      arguments_["combineResultOnly"] = "true";
    }

    if (vm.count("output-raw-prsm-table")) {
      arguments_["outputRawPrsmTable"] = "true";
    }

    if (vm.count("output-prsm-coverage")) {
      arguments_["outputPrsmCoverage"] = "true";
    }

    if (vm.count("disable-post-match")) {
      arguments_["postMassMatch"] = "false";
    }

    // hidden: the PTM localization of -B by the MIScore candidate search of
    // stat/local (the method before toppic_test experiment 140) instead of
    // the default variable PTM search of stat/local_align; --local-align
    // (the old switch for the default) is still accepted
    if (vm.count("local-original")) {
      arguments_["localAlign"] = "false";
    }
    if (vm.count("local-align")) {
      arguments_["localAlign"] = "true";
    }

    if (vm.count("post-min-peak-num")) {
      arguments_["postMinPeakNum"] = post_min_peak_num;
    }
    if (vm.count("average-mass")) {
      arguments_["averageMass"] = "true";
    }
    if (vm.count("disable-zero-shift-prec-error")) {
      arguments_["zeroShiftPrecOneDaError"] = "false";
    }
    if (vm.count("enable-one-shift-prec-error")) {
      arguments_["oneShiftPrecOneDaError"] = "true";
    }
    if (vm.count("enable-var-ptm-prec-error")) {
      arguments_["varPtmPrecOneDaError"] = "true";
    }
    if (vm.count("disable-frag-mass-extension")) {
      arguments_["fragMassExtension"] = "false";
    }
    if (vm.count("frag-mass-dedup")) {
      arguments_["fragMassDedup"] = "true";
    }
    if (vm.count("frag-extend-min-mass")) {
      arguments_["fragExtendMinMass"] = frag_extend_min_mass;
    }
    if (vm.count("var-ptm-filter-cand-num")) {
      arguments_["varPtmFilterCandNum"] = var_ptm_filter_cand_num;
    }
    // hidden: with a deamidation in the variable PTM file, choose between a
    // spectrum's zero-PTM and deamidation PrSMs of the same protein and
    // positions before the top PrSM selection (prsm_deamid_selector,
    // toppic_test experiment 123)
    if (vm.count("deamid-selector")) {
      arguments_["deamidSelector"] = "true";
    }
    // hidden: by default the variable PTM PrSMs' E-values come from a
    // generating function whose residues are the residues and their forms
    // with the variable PTMs, a residue's frequency shared equally by its
    // forms (toppic_test experiment 129); this option gives them the
    // residues without variable PTMs, as every other PrSM
    if (vm.count("disable-var-ptm-evalue-mod-residues")) {
      arguments_["varPtmEvalueModResidues"] = "false";
    }
    if (vm.count("multi-skip-evalue")) {
      arguments_["multiSkipEvalue"] = multi_skip_evalue;
    }
    // hidden: the entries of the filters' mass list
    // (prm_ms_util::getIntMassErrorList), for the one-shift, variable PTM and
    // multi-shift filters, the shift searches' diagonal finder and TopMG's
    // one-shift filter: "merge", one per integer mass (main up to
    // 7b5112c2); "fragment", one per integer mass and base type (the
    // default); "single", one per peak (toppic_test experiment 142)
    if (vm.count("filter-mass-list")) {
      arguments_["filterMassList"] = filter_mass_list;
    }
    // hidden: the tag filter (filter/tag), whose candidate proteins join
    // the lists of the named stages (toppic_test experiment 150)
    if (vm.count("tag-filter")) {
      arguments_["tagFilterStages"] = tag_filter_stages;
    }
    if (vm.count("tag-filter-length")) {
      arguments_["tagFilterLength"] = tag_filter_length;
    }
    if (vm.count("tag-filter-cand-num")) {
      arguments_["tagFilterCandNum"] = tag_filter_cand_num;
    }

    if (vm.count("filtering-result-number")) {
      arguments_["filteringResultNumber"] = filtering_result_num;
    }
  } catch (std::exception& e) {
    std::cerr << "Unhandled Exception in parsing command line" << e.what()
              << ", application will now exit" << std::endl;
    return false;
  }

  return validateArguments();
}

bool ToppicArgument::validateArguments() {
  if (!std::filesystem::exists(arguments_["resourceDir"])) {
    LOG_ERROR("Resource directory "
              << arguments_["resourceDir"] << " does not exist!\n"
              << "Please check if the file directory or name contains special "
                 "characters such as spaces or quotation marks.");
    return false;
  }

  if (!std::filesystem::exists(arguments_["oriDatabaseFileName"])) {
    LOG_ERROR("Database file "
              << arguments_["oriDatabaseFileName"] << " does not exist!\n"
              << "Please check if file directory or name contains special "
                 "characters such as spaces or quotation marks.");
    return false;
  }

  if (!str_util::endsWith(arguments_["oriDatabaseFileName"], ".fasta") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".fa") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".FASTA") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".FA")) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"]
                               << " is not a fasta file!");
    return false;
  }

  if (arguments_["oriDatabaseFileName"].length() > 200) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"]
                               << " path is too long!");
    return false;
  }

  // A database whose protein ids (the header text up to the first
  // whitespace) repeat cannot be indexed by htslib: search its
  // _dedup.fasta copy, which keeps the first record of every id, instead.
  arguments_["oriDatabaseFileName"] =
      fasta_util::removeDuplicatedIds(arguments_["oriDatabaseFileName"]);

  if (arguments_["searchType"] == "TARGET+DECOY") {
    arguments_["databaseFileName"] =
        arguments_["oriDatabaseFileName"] + "_target_decoy";
  } else {
    arguments_["databaseFileName"] =
        arguments_["oriDatabaseFileName"] + "_target";
  }

  for (size_t k = 0; k < spec_file_list_.size(); k++) {
    if (!std::filesystem::exists(spec_file_list_[k])) {
      LOG_ERROR(spec_file_list_[k]
                << " does not exist!\nPlease check if file directory contains "
                   "unproper characters such as spaces/quotation marks");
      return false;
    }

    if (!str_util::endsWith(spec_file_list_[k], ".msalign")) {
      LOG_ERROR("Spectrum file " << spec_file_list_[k]
                                 << " is not a msalign file!");
      return false;
    }

    // Post mass matching (the default) reads the centroided MS/MS peaks from
    // the TopFD SQLite database; fail here rather than after the search.
    if (arguments_["postMassMatch"] == "true") {
      std::string sql_file_name =
          prsm_post_mass_match::sqlFileName(spec_file_list_[k]);
      if (!std::filesystem::exists(sql_file_name)) {
        LOG_ERROR("The TopFD SQLite database "
                  << sql_file_name << " for " << spec_file_list_[k]
                  << " is missing! Post mass matching needs the centroided "
                     "MS/MS peaks stored in it: run TopFD without disabling "
                     "its SQLite database output (without -N) to deconvolute "
                     "the spectra, or run TopPIC with --disable-post-match.");
        return false;
      }
    }

    if (spec_file_list_[k].length() > 200) {
      LOG_ERROR("Spectrum file " << spec_file_list_[k] << " path is too long!");
      return false;
    }

    if (str_util::endsWith(spec_file_list_[k], "_ms1.msalign")) {
      std::cerr << "Warning: Please make sure " << spec_file_list_[k]
                << " is the ms2 spectral file." << std::endl;
    }
  }

  std::string activation = arguments_["activation"];
  if (activation != "CID" && activation != "HCD" && activation != "MPD" &&
      activation != "ETD" && activation != "ETCID" && activation != "ETHCD" &&
      activation != "FILE" && activation != "UVPD") {
    LOG_ERROR("Activation type "
              << activation
              << " error! The value should be "
                 "CID|ETD|ETCID|ETHCD|HCD|MPD|UVPD|FILE!");
    return false;
  }

  std::string search_type = arguments_["searchType"];
  if (search_type != "TARGET" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Search type "
              << search_type
              << " error! The value should be TARGET|TARGET+DECOY!");
    return false;
  }

  // the rescoring trains on the decoys and replaces the E-value in the FDR:
  // off without a target-decoy search (-d) and with --fdr-groups, whose
  // separate FDRs rank by the E-value
  if (search_type != "TARGET+DECOY" || arguments_["useFdrGroups"] == "true") {
    arguments_["rescore"] = "false";
  }

  if (arguments_["secondPrecursor"] == "true") {
    // the second pass runs on one spectrum per set and its PrSMs are
    // merged with the first pass's before the report of the file
    if (arguments_["combinedOutputName"] != "") {
      LOG_ERROR("--second-precursor cannot be combined with -c!");
      return false;
    }
    if (arguments_["groupSpectrumNumber"] != "1") {
      LOG_ERROR("--second-precursor needs one spectrum per set (-r 1)!");
      return false;
    }
  }

  std::string shift_number = arguments_["shiftNumber"];
  if (shift_number != "0" && shift_number != "1" && shift_number != "2") {
    LOG_ERROR("Unexpected modification number "
              << shift_number << " error! The value should be 0|1|2!");
    return false;
  }

  std::string allow_mod = arguments_["allowProtMod"];
  std::vector<std::string> strs = str_util::split(allow_mod, ",");
  for (size_t i = 0; i < strs.size(); i++) {
    if (strs[i] != "NONE" && strs[i] != "M_ACETYLATION" && strs[i] != "NME" &&
        strs[i] != "NME_ACETYLATION") {
      LOG_ERROR(
          "N-Terminal Variable modifications can only be NONE, M_ACETYLATION, "
          "NME or NME_ACETYLATION.");
      return false;
    }
  }

  std::string variable_ptm_num = arguments_["variablePtmNum"];
  try {
    int num = std::stoi(variable_ptm_num.c_str());
    if (num < 0) {
      LOG_ERROR("Maximum variable modification number "
                << num << " error! The value should be nonnegative.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Maximum variable modification number "
              << variable_ptm_num << " should be a number.");
    return false;
  }

  std::string cutoff_spectral_type = arguments_["cutoffSpectralType"];
  if (cutoff_spectral_type != "EVALUE" && cutoff_spectral_type != "FDR") {
    LOG_ERROR("Spectrum-level cutoff type "
              << cutoff_spectral_type
              << " error! The value should be EVALUE|FDR");
    return false;
  }

  if (cutoff_spectral_type == "FDR" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Spectrum-level cutoff type "
              << cutoff_spectral_type
              << " error! FDR cutoff cannot be used when no decoy database is "
                 "used! Please add argument '-d' in the command.");
    return false;
  }

  std::string cutoff_proteoform_type = arguments_["cutoffProteoformType"];
  if (cutoff_proteoform_type != "EVALUE" && cutoff_proteoform_type != "FDR") {
    LOG_ERROR("Proteoform-level cutoff type "
              << cutoff_proteoform_type
              << " error! The value should be EVALUE|FDR");
    return false;
  }

  if (cutoff_proteoform_type == "FDR" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Proteoform-level cutoff type "
              << cutoff_proteoform_type
              << " error! FDR cutoff cannot be used when no decoy database is "
                 "used! Please add argument '-d' in the command.");
    return false;
  }

  std::string cutoff_protein_type = arguments_["cutoffProteinType"];
  if (cutoff_protein_type != "EVALUE" && cutoff_protein_type != "FDR") {
    LOG_ERROR("Protein-level cutoff type "
              << cutoff_protein_type
              << " error! The value should be EVALUE|FDR");
    return false;
  }

  if (cutoff_protein_type == "FDR" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Protein-level cutoff type "
              << cutoff_protein_type
              << " error! FDR cutoff cannot be used when no decoy database is "
                 "used! Please add argument '-d' in the command.");
    return false;
  }

  std::string min_shift_mass = arguments_["minShiftMass"];
  try {
    std::stod(min_shift_mass);
  } catch (const std::exception& ex) {
    LOG_ERROR("Minimum modification mass " << min_shift_mass
                                           << " should be a number.");
    return false;
  }

  std::string max_shift_mass = arguments_["maxShiftMass"];
  try {
    double mass = std::stod(max_shift_mass);
    if (mass <= 0.0) {
      LOG_ERROR("Maximum modification mass "
                << max_shift_mass << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Maximum modification mass " << max_shift_mass
                                           << " should be a number.");
    return false;
  }

  std::string n_term_label_mass = arguments_["nTermLabelMass"];
  try {
    double mass = std::stod(n_term_label_mass);
    if (mass < 0.0) {
      LOG_ERROR("N-terminal label mass "
                << n_term_label_mass
                << " error! The value should not be negative.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("N-terminal label mass " << n_term_label_mass
                                       << " should be a number.");
    return false;
  }

  for (const auto& [key, name] :
       {std::pair<const char*, const char*>{"precErrorTolerance",
                                            "Precursor error tolerance"},
        std::pair<const char*, const char*>{"fragErrorTolerance",
                                            "Fragment error tolerance"}}) {
    std::string tole_value = arguments_[key];
    try {
      double tole = std::stoi(tole_value);
      if (tole <= 0) {
        LOG_ERROR(name << ": " << tole_value
                       << " error! The value should be positive.");
        return false;
      }
    } catch (const std::exception& ex) {
      LOG_ERROR(name << ": " << tole_value << " should be a number.");
      return false;
    }
  }

  std::string prec_dalton_error = arguments_["precDaltonError"];
  if (prec_dalton_error != "") {
    bool valid = false;
    try {
      std::size_t end = 0;
      double tole = std::stod(prec_dalton_error, &end);
      valid = end == prec_dalton_error.size() && std::isfinite(tole) &&
              tole > 0;
    } catch (const std::exception& ex) {
      valid = false;
    }
    if (!valid) {
      LOG_ERROR("Precursor error tolerance in Dalton (--prec-dalton-error) "
                << prec_dalton_error << " should be a positive number.");
      return false;
    }
  }

  std::string env_cnn_cutoff_value = arguments_["envCnnCutoff"];
  try {
    double cutoff = std::stod(env_cnn_cutoff_value);
    if (cutoff < 0 || cutoff > 1) {
      LOG_ERROR("EnvCNN score cutoff: "
                << env_cnn_cutoff_value
                << " error! The value should be between 0 and 1.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("EnvCNN score cutoff: " << env_cnn_cutoff_value
                                      << " should be a number.");
    return false;
  }

  std::string form_error_tole_value = arguments_["proteoformErrorTolerance"];
  try {
    double tole = std::stod(form_error_tole_value);
    if (tole <= 0) {
      LOG_ERROR("PrSM clustering error tolerance: "
                << form_error_tole_value
                << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("PrSM clustering error tolerance: " << form_error_tole_value
                                                  << " should be a number.");
    return false;
  }

  std::string cutoff_spectral_value = arguments_["cutoffSpectralValue"];
  try {
    double th = std::stod(cutoff_spectral_value);
    if (th < 0) {
      LOG_ERROR("Spectrum-level cutoff value "
                << cutoff_spectral_value
                << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Spectrum-level cutoff value " << cutoff_spectral_value
                                             << " should be a number.");
    return false;
  }

  std::string cutoff_proteoform_value = arguments_["cutoffProteoformValue"];
  try {
    double th = std::stod(cutoff_proteoform_value);
    if (th < 0) {
      LOG_ERROR("Proteoform-level cutoff value "
                << cutoff_proteoform_value
                << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Proteoform-level cutoff value " << cutoff_proteoform_value
                                               << " should be a number.");
    return false;
  }

  std::string cutoff_protein_value = arguments_["cutoffProteinValue"];
  try {
    double th = std::stod(cutoff_protein_value);
    if (th < 0) {
      LOG_ERROR("Protein-level cutoff value "
                << cutoff_protein_value
                << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Protein-level cutoff value " << cutoff_protein_value
                                            << " should be a number.");
    return false;
  }

  std::string thread_number = arguments_["threadNumber"];
  try {
    int num = std::stoi(thread_number);
    bool valid = mem_check::checkThreadNum(num, "toppic");
    if (!valid) {
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Thread number " << thread_number << " should be a number.");
    return false;
  }

  std::string post_min_peak_num = arguments_["postMinPeakNum"];
  try {
    if (std::stoi(post_min_peak_num) < 1) {
      LOG_ERROR("Post mass matching minimum isotopic peak number "
                << post_min_peak_num << " should be at least 1.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Post mass matching minimum isotopic peak number "
              << post_min_peak_num << " should be a number.");
    return false;
  }

  std::string frag_extend_min_mass = arguments_["fragExtendMinMass"];
  try {
    if (std::stod(frag_extend_min_mass) < 0) {
      LOG_ERROR("Fragment extension minimum mass "
                << frag_extend_min_mass << " should not be negative.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Fragment extension minimum mass " << frag_extend_min_mass
                                                 << " should be a number.");
    return false;
  }

  // the variable PTM filter's candidates kept per spectrum, for complete,
  // prefix and suffix, and internal proteoforms (VarPtmFilterMng)
  std::string var_cand_num = arguments_["varPtmFilterCandNum"];
  if (var_cand_num != "") {
    std::vector<std::string> nums = str_util::split(var_cand_num, ",");
    bool valid = nums.size() == 3;
    for (size_t i = 0; valid && i < nums.size(); i++) {
      valid = !nums[i].empty() &&
              nums[i].find_first_not_of("0123456789") == std::string::npos &&
              nums[i].size() <= 6 && std::stoi(nums[i]) > 0;
    }
    if (!valid) {
      LOG_ERROR("--var-ptm-filter-cand-num "
                << var_cand_num
                << " should be three positive integers separated by commas "
                   "(complete,prefix/suffix,internal).");
      return false;
    }
  }

  if (arguments_["evalueComplementTolerance"] != "strict" &&
      arguments_["evalueComplementTolerance"] != "relax" &&
      arguments_["evalueComplementTolerance"] != "separate") {
    LOG_ERROR("--evalue-complement-tolerance "
              << arguments_["evalueComplementTolerance"]
              << " is not strict, relax or separate");
    return false;
  }
  if (!SpPara::isFilterMassListModeName(arguments_["filterMassList"])) {
    LOG_ERROR("--filter-mass-list " << arguments_["filterMassList"]
                                    << " error! The value should be merge, "
                                       "fragment or single!");
    return false;
  }

  // the tag filter: the stages, the tag length (at least the index's K-mer
  // length) and the candidates per spectrum
  std::set<TagFilterMng::Stage> tag_stages;
  if (!TagFilterMng::parseStages(arguments_["tagFilterStages"], &tag_stages)) {
    LOG_ERROR("--tag-filter " << arguments_["tagFilterStages"]
                              << " error! The value should be a "
                                 "comma-separated list of var, one, multi "
                                 "and topmg, or all!");
    return false;
  }
  auto isPositiveInt = [](const std::string& str, int min_value) {
    return !str.empty() && str.size() <= 6 &&
           str.find_first_not_of("0123456789") == std::string::npos &&
           std::stoi(str) >= min_value;
  };
  if (!isPositiveInt(arguments_["tagFilterLength"], TagIndex::K)) {
    LOG_ERROR("--tag-filter-length " << arguments_["tagFilterLength"]
                                     << " error! The value should be an "
                                        "integer of at least "
                                     << TagIndex::K << ".");
    return false;
  }
  if (!isPositiveInt(arguments_["tagFilterCandNum"], 1)) {
    LOG_ERROR("--tag-filter-cand-num " << arguments_["tagFilterCandNum"]
                                       << " error! The value should be a "
                                          "positive integer.");
    return false;
  }

  std::string topmg_mode = arguments_["topmgMode"];
  if (topmg_mode != "disabled" && topmg_mode != "pipeline" &&
      topmg_mode != "topmg-only") {
    LOG_ERROR("--topmg " << topmg_mode
                         << " error! The value should be disabled, pipeline "
                            "or topmg-only!");
    return false;
  }
  // TopMG's identification (hidden --topmg pipeline or topmg-only) searches
  // whole proteins or whole proteins and their fragments; no other
  // combination of the proteoform types is meaningful for its proteoform
  // graphs
  if (topmg_mode != "disabled") {
    std::string prot_type = arguments_["allowProtType"];
    if (prot_type != "COMPLETE,PREFIX,SUFFIX,INTERNAL" &&
        prot_type != "COMPLETE") {
      LOG_ERROR("Proteoform type "
                << prot_type
                << " error! With TopMG's identification, the value should be "
                   "COMPLETE,PREFIX,SUFFIX,INTERNAL or COMPLETE!");
      return false;
    }
  }
  // TopMG's identification needs the variable modifications; as the only
  // search, a run without them would search nothing
  if (topmg_mode == "topmg-only" && arguments_["variablePtmFileName"] == "") {
    LOG_ERROR("--topmg topmg-only error! It needs a variable modification "
              "file (-b).");
    return false;
  }
  std::string topmg_align_error_tole = arguments_["topmgAlignErrorTolerance"];
  try {
    if (std::stod(topmg_align_error_tole) <= 0) {
      LOG_ERROR("--topmg-align-error-tolerance " << topmg_align_error_tole
                << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("--topmg-align-error-tolerance " << topmg_align_error_tole
              << " should be a number.");
    return false;
  }

  // the E-value method of TopMG's identification means nothing without it
  if (arguments_["topmgUseMcmc"] == "true" && topmg_mode == "disabled") {
    LOG_ERROR("--topmg-use-mcmc error! It needs --topmg pipeline or "
              "topmg-only.");
    return false;
  }

  // empty: the default of the multiple-shift stage (the smallest decoy
  // E-value with a decoy database, 0.01 without)
  std::string multi_skip_evalue = arguments_["multiSkipEvalue"];
  if (multi_skip_evalue != "") {
    try {
      if (std::stod(multi_skip_evalue) < 0) {
        LOG_ERROR("Multiple-shift skip E-value " << multi_skip_evalue
                                                 << " should not be negative.");
        return false;
      }
    } catch (const std::exception& ex) {
      LOG_ERROR("Multiple-shift skip E-value " << multi_skip_evalue
                                               << " should be a number.");
      return false;
    }
  }

  return true;
}

}  // namespace toppic
