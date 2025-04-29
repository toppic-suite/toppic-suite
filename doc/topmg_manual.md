## TopMG
TopMG is a software tool for identifying highly modified proteoforms by searching top-down tandem mass spectra against a protein sequence database. It is capable of identifying proteoforms with multiple variable PTMs and unexpected alterations, such as histone proteoforms and phosphorylated ones. It uses mass graphs, which efficiently represent candidate proteoforms with multiple variable PTMs, to increase the speed and sensitivity in proteoform identification. In addition, approximate spectrum-based filtering methods are employed for protein sequence filtering, and a Markov chain Monte Carlo method (TopMCMC) is used for estimating the statistical significance of identifications.

### Minimum requirements
A computer with at least 16 GB memory and a 64-bit Linux or Windows operating system is required. 

### Input

* A protein database file in the FASTA format
* A mass spectrum data file in the msalign format
* A text file of variable PTMs
* A text file of fixed PTMs (optional)
* A text file containing LC-MS feature information (optional)

### Output

TopMG outputs two TSV files, an XML file, and a collection of HTML files for identified proteoforms. For example, when the input mass spectrum data file is spectra_ms2.msalign, the output includes:

* **spectra_ms2_topmg_prsm.tsv**: a TSV file containing identified PrSMs with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, all the proteins are reported.

* **spectra_ms2_topmg_prsm_single.tsv**: a TSV file containing identified PrSMs with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.

* **spectra_ms2_topmg_proteoform.tsv**: a TSV file containing identified proteoforms with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, all the proteins are reported.

* **spectra_ms2_topmg_proteoform_single.tsv**: a TSV file containing identified proteoforms with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.

* **spectra_ms2_topmg_proteoform.xml**: an XML file containing identified proteoforms with the E-value or proteoform-level FDR cutoff.

* **spectra_ms2_topmg_prsm.xml**: an XML file containing all identified PrSMs without clustering and filtering.

* **spectra_html/topmg_prsm_cutoff**: a folder containing JavaScript files of identified PrSMs using the E-value or spectrum-level FDR cutoff.

* **spectra_html/topmg_proteoform_cutoff**: a folder containing JavaScript files of identified PrSMs using the E-value or proteoform-level cutoff.

* **spectra_html/topmsv**: a folder containing HTML files for the visualization of identified PrSMs.

To browse identified proteins, proteoforms, and PrSMs, use a chrome browser to open the file spectra_html/topmsv/index.html. Google Chrome is recommended (Firefox and Edge are not recommended).

When the input contains two or more spectrum files, TopMG outputs two TSV files, an XML file, and a collection of HTML files for each input file. When a file name is specified for combined identifications, it combines spectra and proteoforms identified from all the input files, removes redundant proteoform identifications, and reports two TSV files, an XML file, and a collection of HTML files for the combined results. For example, when the input is spectra1_ms2.msalign and spectra2_ms2.msalign and the combined output file name is "combined," the output files are:

* **combined_ms2_topmg_prsm.tsv**: a TSV file containing PrSMs identified from all the input files with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, all the proteins are reported.

* **combined_ms2_topmg_prsm_single.tsv**: a TSV file containing PrSMs identified from all the input files with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.

* **combined_ms2_topmg_proteoform.tsv**: a TSV file containing proteoforms identified from all the input files with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, all the proteins are reported.

* **combined_ms2_topmg_proteoform_single.tsv**: a TSV file containing proteoforms identified from all the input files with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.

* **combined_ms2_topmg_proteoform.xml**: an XML file containing proteoforms identified from all the input files with the E-value or proteoform-level FDR cutoff.

* **combined_ms2_topmg_prsm.xml**: an XML file containing all identified PrSMs without clustering and filtering.

* **combined_html/topmg_prsm_cutoff**: a folder containing JavaScript files of PrSMs identified from all the input files using the E-value or spectrum-level FDR cutoff.

* **combined_html/topmg_proteoform_cutoff**: a folder containing JavaScript files of PrSMs identified from all the input files using the E-value or proteoform-level cutoff.

* **combined_html/topmsv**: a folder containing HTML files for the visualization of identified PrSMs.

### Command line usage

To run TopMG, open a terminal window and run the following command.
```
topmg [options] database-file-name spectrum-file-names
```

Options
```
-h [ --help ]
```
Print the help message.
```
-a [ --activation ] <CID|HCD|ETD|UVPD|FILE>
```
Fragmentation method of MS/MS spectra. When FILE is used, fragmentation methods of spectra are given in the input spectral data file. Default value: FILE.
```
-f [ --fixed-mod ] <C57|C58|a fixed modification file>
```
Set fixed modifications. Three available options: C57, C58, or the name of a text file specifying fixed modifications (see an example file). When C57 is selected, carbamidomethylation on cysteine is the only fixed modification. When C58 is selected, carboxymethylation on cysteine is the only fixed modification.
```
-n [ --n-terminal-form ] <a list of allowed N-terminal forms>
```
Set N-terminal forms of proteins. Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, and M_ACETYLATION. NONE stands for no modifications, NME for N-terminal methionine excision, NME_ACETYLATION for N-terminal acetylation after the initiator methionine is removed, and M_ACETYLATION for N-terminal methionine acetylation. When multiple forms are allowed, they are separated by commas. Default value: NONE,M_ACETYLATION,NME,NME_ACETYLATION.
```
-d [ --decoy ]
```
Use a shuffled decoy protein database to estimate spectrum and proteoform level FDRs. When -d is chosen, a shuffled decoy database is automatically generated and appended to the target database before database search, and FDR rates are estimated using the target-decoy approach.
```
-e [ --mass-error-tolerance ] <a positive integer>
```
Set the error tolerance for precursor and fragment masses in ppm. Default value: 10 ppm.
```
-p [ --proteoform-error-tolerance ] <a positive number>
```
Set the error tolerance for identifying PrSM clusters (in Dalton). Default value: 1.2 Dalton.
```
-M [ --max-shift ] <a number>
```
Set the maximum absolute value for unexpected mass shifts (in Dalton). Default value: 500 Dalton.
```
-t [ --spectrum-cutoff-type ] <EVALUE|FDR>
```
Set the spectrum level cutoff type for filtering PrSMs. Default value: EVALUE.
```
-v [ --spectrum-cutoff-value ] <a positive number>
```
Set the spectrum level cutoff value for filtering PrSMs. Default value: 0.01.
```
-T [ --proteoform-cutoff-type ] <EVALUE|FDR>
```
Set the proteoform level cutoff type for filtering proteoforms and PrSMs. Default value: EVALUE.
```
-V [ --proteoform-cutoff-value ] <a positive number>
```
Set the proteoform level cutoff value for filtering proteoforms and PrSMs. Default value: 0.01.
```
-i [ --mod-file-name ] <a modification file>
```
Specify a text file of variable PTMs. See an example file.
```
-u [ --thread-number ] <a positive number>
```
Set the number of threads used in the computation. Default value: 1.
```
-x [ --no-topfd-feature ]
```
Specify that there are no TopFD feature files for proteoform identification.
```
-D [ --use-asf-diagonal ]
```
Use the ASF-DIAGONAL method for protein sequence filtering. The default filtering method is ASF-RESTRICT. When -D is selected, both ASF-RESTRICT and ASF-DIAGONAL will be used. The combined approach may identify more PrSMs, but it is much slower than using ASF-RESTRICT only. See this paper for more details.
```
-P [ --var-ptm ] <a positive number>
```
Set the maximum number of variable PTM sites in a proteoform. Default value: 5.
```
-s --num-shift <0|1|2>
```
Set the maximum number of unexpected mass shifts in a proteoform. Default value: 0.
```
-w [ --whole-protein-only ]
```
Report only proteoforms of whole protein sequences.
```
-c [ --combined-file-name ] <a filename>
```
Specify an output file name for combined identifications when the input consists of multiple spectrum files.
```
-k [ --keep ]
```
Keep intermediate files.
```
-K [ --keep-decoy-ids ]
```
Keep decoy identifications.
```
-g [ --skip-html-folder ]
```
Skip the generation of HTML files for visualization.
Advanced options
```
-j [ --proteo-graph-dis ] <a positive number>
```
Set the length of the largest gap in constructing proteoform graphs. Default value: 40. See this paper for more details.
```
-G [ --var-ptm-in-gap ] <a positive number>
```
Set the maximum number of variable PTM sites in a gap in a proteoform graph. Default value: 5. See this paper for more details.

### Examples

To use the following examples, the current folder needs to contain a variable modification file variable_mods.txt. (See an example.)

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file spectra_ms2.feature. The user does not need to specify the feature file name. TopMG will automatically obtain the feature file name from the spectrum file name spectra_ms2.msalign.
```
topmg -i variable_mods.txt proteins.fasta spectra_ms2.msalign
```
Search two deconvoluted MS/MS spectrum files spectra1_ms2.msalign and spectra2_ms2.msalign against a protein database file proteins.fasta with feature files. In addition, all identifications are combined and reported using a file name "combined."
```
topmg -i variable_mods.txt -c combined proteins.fasta spectra1_ms2.msalign spectra2_ms2.msalign
```
Search all deconvoluted MS/MS spectrum files in the current folder against a protein database file proteins.fasta with feature files.
```
topmg -i variable_mods.txt proteins.fasta *_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta without feature files.
```
topmg -i variable_mods.txt -x proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file and a fixed modification: carbamidomethylation on cysteine.
```
topmg -i variable_mods.txt -f C57 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. In an identified proteoform, at most 1 unexpected mass shift and 4 variable PTMs are allowed and the maximum value for unexpected mass shifts is 10,000 Dalton.
```
topmg -i variable_mods.txt -P 4 -s 1 -M 10000 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. The error tolerance for precursor and fragment masses is 5 ppm.
```
topmg -i variable_mods.txt -e 5 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. Use the target decoy approach to compute spectrum level and proteoform level FDRs, filter identified proteoform spectrum-matches by a 5% spectrum-level FDR, and filter identified proteoforms by a 5% proteoform-level FDR.
```
topmg -i variable_mods.txt -d -t FDR -v 0.05 -T FDR -V 0.05 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. Use 6 CPU threads to speed up the computation.
```
topmg -i variable_mods.txt -u 6 proteins.fasta spectra_ms2.msalign
```