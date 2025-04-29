## TopPIC

TopPIC identifies and characterizes proteoforms at the proteome level by searching top-down tandem mass spectra against a protein sequence database. It efficiently identifies proteoforms with post-translational modificatons (PTMs) and unexpected alterations, such as mutations, accurately estimates the statistical significance of identifications, and characterizes reported proteoforms with unknown mass shifts. It uses several techniques, such as indexes, spectral alignment, generating function methods, and the modification identification score (MIScore), to increase the speed, sensitivity, and accuracy.

### Minimum requirements
A computer with at least 16 GB memory and a 64-bit Linux or Windows operating system is required. 

### Input

* A protein database file in the FASTA format
* A mass spectrum data file in the msalign format
* A text file containing LC-MS feature information (optional)
* A text file of fixed PTMs (optional)
* A text file of variable PTMs (optional)
* A text file of PTMs for the characterization of unexpected mass shifts (optional)

### Output

TopPIC outputs four tab separated value (TSV) files, two XML files, and a collection of HTML files for identified proteoforms. For example, when the input data file is spectra_ms2.msalign, the output includes:

* **spectra_ms2_toppic_prsm.tsv**: a TSV file containing identified proteoform spectrum-matches (PrSMs) with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, all the proteins are reported.

* **spectra_ms2_toppic_prsm_single.tsv**: a TSV file containing identified proteoform spectrum-matches (PrSMs) with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.

* **spectra_ms2_toppic_proteoform.tsv**: a TSV file containing identified proteoforms with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, all the proteins are reported.

* **spectra_ms2_toppic_proteoform_single.tsv**: a TSV file containing identified proteoforms with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.

* **spectra_ms2_toppic_proteoform.xml**: an XML file containing identified proteoforms with the E-value or proteoform-level FDR cutoff.

* **spectra_ms2_toppic_prsm.xml**: an XML file containing all identified PrSMs without clustering and filtering.

* **spectra_html/toppic_prsm_cutoff**: a folder containing JavaScript files of identified PrSMs using the E-value or spectrum-level FDR cutoff.

* **spectra_html/toppic_proteoform_cutoff**: a folder containing JavaScript files of identified PrSMs using the E-value or proteoform-level FDR cutoff.

* **spectra_html/topmsv**: a folder containing HTML files for the visualization of identified PrSMs.

To browse identified proteins, proteoforms, and PrSMs, use a chrome browser to open the file spectrum_html/topmsv/index.html. Google Chrome is recommended (Firefox and Edge are not recommended).

When the input contains two or more data files, TopPIC outputs four TSV files, two XML files, and a collection of HTML files for each input file. When a file name is specified for combined identifications, it combines spectra and proteoforms identified from all the input files, removes redundant proteoform identifications, and reports four TSV files, two XML files, and a collection of HTML files for the combined results. For example, when the input is spectra1_ms2.msalign and spectra2_ms2.msalign and the combined output file name is "combined," the output files are:

* **combined_ms2_toppic_prsm.tsv**: a TSV file containing PrSMs identified from all the input files with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, all the proteins are reported.

* **combined_ms2_toppic_prsm_single.tsv**: a TSV file containing PrSMs identified from all the input files with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.

* **combined_ms2_toppic_proteoform.tsv**: a TSV file containing proteoforms identified from all the input files with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, all the proteins are reported.

* **combined_ms2_toppic_proteoform_single.tsv**: a TSV file containing proteoforms identified from all the input files with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.

* **combined_ms2_toppic_proteoform.xml**: an XML file containing proteoforms identified from all the input files with the E-value or proteoform-level FDR cutoff.

* **combined_ms2_toppic_prsm.xml**: an XML file containing all identified PrSMs without clustering and filtering.

* **combined_html/toppic_prsm_cutoff**: a folder containing JavaScript files of PrSMs identified from all the input files using the E-value or spectrum-level FDR cutoff.

* **combined_html/toppic_proteoform_cutoff**: a folder containing JavaScript files of PrSMs identified from all the input files using the E-value or proteoform-level FDR cutoff.

* **combined_html/topmsv**: a folder containing HTML files for the visualization of identified PrSMs.

### Command line usage

To run TopPIC, open a terminal window and run the following command.
```
toppic [options] database-file-name spectrum-file-names
```

Options
```
-h [ --help ]
```
Print the help message.
```
-a [ --activation ] <CID|HCD|ETD|UVPD|FILE>
```
Set the fragmentation method(s) of MS/MS spectra. When "FILE" is selected, the fragmentation methods of spectra are given in the input spectrum data file. Default value: FILE.
```
-f [ --fixed-mod ] <C57|C58|a fixed modification file>
```
Set fixed modifications. Three available options: C57, C58, or the name of a text file containing the information of fixed modifications (see an example file). When C57 is selected, carbamidomethylation on cysteine is the only fixed modification. When C58 is selected, carboxymethylation on cysteine is the only fixed modification.
```
-n [ --n-terminal-form ] <a list of allowed N-terminal forms>
```
Set N-terminal forms of proteins. Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, and M_ACETYLATION. NONE stands for no modifications, NME for N-terminal methionine excision, NME_ACETYLATION for N-terminal acetylation after the initiator methionine is removed, and M_ACETYLATION for N-terminal methionine acetylation. When multiple forms are allowed, they are separated by commas. Default value: NONE,M_ACETYLATION,NME,NME_ACETYLATION.
```
-s [ --num-shift ] <0|1|2>
```
The maximum number of unexpected mass shifts in a PrSM. Default value: 1.
```
-m [ --min-shift ] <a number>
```
The minimum value for unexpected mass shifts (in Dalton). Default value: -500 Dalton.
```
-M [ --max-shift ] <a number>
```
The maximum value for unexpected mass shifts (in Dalton). Default value: 500 Dalton.
```
-S [ --variable-ptm-num] <a number>
```
The maximum number of variable PTM sites in a proteoform-spectrum-match. Default value: 3.
```
-b [ --variable-ptm-file-name] a variable PTM file
```
Specify a text file containing the information of varaible PTMs (see an example variable PTM file).
```
-d [ --decoy ]
```
Use a shuffled decoy protein database to estimate spectrum and proteoform-level FDRs. When -d is chosen, a shuffled decoy database is automatically generated and appended to the target database before database search, and FDRs are estimated using the target-decoy approach.
```
-e [ --mass-error-tolerance ] <a positive integer>
```
Set the error tolerance for precursor and fragment masses in part-per-million (ppm). Default value: 10.
```
-p [ --proteoform-error-tolerance ] <a positive number>
```
Set the error tolerance for identifying PrSM clusters (in Dalton). Default value: 1.2 Dalton.
```
-t [ --spectrum-cutoff-type ] <EVALUE|FDR>
```
Set the spectrum-level cutoff type for filtering PrSMs. Default value: EVALUE.
```
-v [ --spectrum-cutoff-value ] <a positive number>
```
Set the spectrum-level cutoff value for filtering PrSMs. Default value: 0.01.
```
-T [ --proteoform-cutoff-type ] <EVALUE|FDR>
```
Set the proteoform-level cutoff type for filtering proteoforms and PrSMs. Default value: EVALUE.
```
-V [ --proteoform-cutoff-value ] <a positive number>
```
Set the proteoform-level cutoff value for filtering proteoforms and PrSMs. Default value: 0.01.
```
-A [ --approximate-spectra ]
```
Use approximate spectra to increase the sensitivity in protein filtering (see this paper for details).
```
-l [ --lookup-table ]
```
Use a lookup table method for computing E-values. It is faster than the default generating function approach, but it may reduce the number of identifications.
```
-B [ --local-ptm-file-name ] <a common modification file>
```
Specify a text file containing a list of common PTMs for proteoform characterization. The PTMs are used to identify and localize PTMs in reported PrSMs with unknown mass shifts. See an example file.
```
-H [ --miscore-threshold ] <a number between 0 and 1>
```
Set the MIScore threshold (see paper) for filtering results of PTM characterization. Default value: 0.15.
```
-u [ --thread-number ] <a positive number>
```
Set the number of threads used in the computation. Default value: 1.
```
-r [ --num-combined-spectra ] <a positive integer>
```
Set the number of combined spectra. The parameter is set to 2 (or 3) for combining spectral pairs (or triplets) generated by the alternating fragmentation mode. Default value: 1.
```
-c [ --combined-file-name ] <a filename>
```
Specify an output file name for combined identifications when the input consists of multiple data files.
```
-x [ --no-topfd-feature ]
```
Specify that there are no TopFD feature files for proteoform identification.
```
-k [ --keep-temp-files ]
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

### Examples

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file spectra_ms2.feature (reported by TopFD). The user does not need to specify the feature file name. TopPIC will automatically obtain the feature file name from the spectrum file name spectra_ms2.msalign.
```
toppic proteins.fasta spectra_ms2.msalign
```
Search two deconvoluted MS/MS spectrum files spectra1_ms2.msalign and spectra2_ms2.msalign against a protein database file proteins.fasta with feature files. In addition, all identifications are combined and reported using a file name "combined."
```
toppic -c combined proteins.fasta spectra1_ms2.msalign spectra2_ms2.msalign
```
Search all deconvoluted MS/MS spectrum files in the current folder against a protein database file proteins.fasta with feature files.
```
toppic proteins.fasta *_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta without feature files.
```
toppic -x proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file and a fixed modification: carbamidomethylation on cysteine.
```
toppic -f C57 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. In an identified proteoform, at most 2 mass shifts are allowed and the maximum allowed mass shift value is 10,000 Dalton.
```
toppic -s 2 -M 10000 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. Two variable PTMs: oxidation on M and methylation on K are used. The modification file two_var_mods.txt can be found here.
```
toppic -b two_var_mods.txt proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. The error tolerance for precursor and fragment masses is 5 ppm.
```
toppic -e 5 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. Use the target decoy approach to compute spectrum level and proteoform-level FDRs, filter identified proteoform spectrum-matches by a 5% spectrum level FDR, and filter identified proteoforms by a 5% proteoform-level FDR.
```
toppic -d -t FDR -v 0.05 -T FDR -V 0.05 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign with alternating CID, HCD, and ETD spectra against a protein database file proteins.fasta with a feature file. Combine alternating CID, HCD, and ETD spectra to increase proteoform coverage.
```
toppic -r 3 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. After proteoforms with unexpected mass shifts are identified, TopPIC matches the mass shifts to four common PTMs: acetylation, phosphorylation, oxidation and methylation, and uses an MIScore cutoff 0.1 to filter reported PTM sites. The modification file common_mods.txt can be found here.
```
toppic -B common_mods.txt -H 0.1 proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file.Use 6 CPU threads to speed up the computation.
```
toppic -u 6 proteins.fasta spectra_ms2.msalign 
```