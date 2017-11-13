To run TopPIC suite, a computer with at least 4 GB memory and a 64-bit Linux or Windows operating system is required. TopFD and TopPIC provide a command line interface and a graphical user interface (GUI) for both Linux and Windows users.

# TopFD

## Input

The input of TopFD is an mzML or mzXML top-down mass spectrometry data file. Raw mass spectral data generated from various mass spectrometers can be converted to mzML or mzXML files using msconvert.

```
msconvert.exe --filter "peakPicking true 1-" --mzML --zlib spectra.raw
```

## Output

TopFD outputs an MS1 feature text file with a file extension `feature` and two deconvoluted mass spectral data files (one is for MS1 spectra and the other for MS/MS spectra) in the msalign format with a file extension `msalign`, which is similar to the MGF file format. For example, when the input file name is spectra.mzML, the output files include:

* __spectra.feature__: a feature file containing MS/MS scan ids and their corresponding LC/MS features.
* __spectra_ms1.msalign__: a list of deconvoluted MS1 spectra.
* __spectra_ms2.msalign__: a list of deconvoluted MS/MS spectra.

## Command line usage

To run TopFD, open a console and run the following command.
```
topfd [options] spectrum-file-name
```

### Options

```
-h [ --help ]
```
Print the help message.
```
-c [ --max-charge ] <a positive integer>
```
Set the maximum charge state of precursor and fragment ions. The default value is 30.
```
-m [ --max-mass ] <a positive number>
```
Set the maximum monoisotopic mass of precursor and fragment ions. The default value is 100,000 Dalton.
```
-e [ --mz-error ] <a positive number>
```
Set the error tolerance of m/z values of spectral peaks. The default value is 0.02 m/z.
```
-s [ --sn-ratio ] <a positive number>
```
Set the signal/noise ratio. The default value is 1.
```
-w [ --precursor-window ] <a positive number>
```
Set the precursor isolation window size. The default value is 3.0 m/z.
```
-n [ --missing-level-one ]
```
The input spectrum file does not contain MS1 spectra.

# TopPIC

## Input

* A protein database file in the FASTA format
* A mass spectrum data file in the msalign format
* A text file to specify fixed PTMs (optional)
* A text file containing MS1 feature information (optional)
* A text file to specify PTMs for the characterization of mass shifts (optional)

## Output

TopPIC outputs two tab delimited text files and a collection of html files for identified proteoforms. For example, when the input mass spectrum data file is spectra_ms2.msalign, the output includes:

* spectra_ms2.OUTPUT_TABLE: a tab delimited text file containing a list of identified proteoform spectrum-matches at the spectrum level.
* spectra_ms2.FORM_OUTPUT_TABLE: a tab delimited text file containing a list of identified proteoform spectrum-matches at the proteoform level.
* spectra_ms2_prsm_cutoff_html: a folder containing html files for the annotation of identified proteoform spectrum-matches using the spectrum level E-value or FDR cutoff.
* spectra_ms2_proteoform_cutoff_html: a folder containing html files for the annotation of identified proteoform spectrum-matches using the proteoform level E-value or FDR cutoff.

## Command line usage

To run TopPIC, open a console and run the following command.
```
toppic [options] database-file-name spectrum-file-name
```
Options
```
-h [ --help ]
```
Print the help message.
```
-a [ --activation ] <CID|HCD|ETD|UVPD|FILE>
```
Set the fragmentation method(s) of tandem mass spectra. When FILE is used, fragmentation methods of spectra are given in the input spectrum data file. Default value: FILE.
```
-f [ --fixed-mod ] <C57|C58|a fixed modification file>
```
Set the fixed modifications. Three available options: C57, C58, or the name of a text file containing the information of fixed modifications (see example file). When C57 is selected, carbamidomethylation on cysteine is the only fixed modification. When C58 is selected, carboxymethylation on cysteine is the only fixed modification.
```
-n [ --n-terminal-form ] <a list of allowed N-terminal forms>
```
Set the N-terminal forms of proteins. Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, and M_ACETYLATION. NONE stands for no modifications, NME for N-terminal methionine excision, NME_ACETYLATION for N-terminal acetylation after the initiator methionine is removed, and M_ACETYLATION for N-terminal methionine acetylation. When multiple forms are allowed, they are separated by commas. Default value: NONE,M_ACETYLATION,NME,NME_ACETYLATION.
```
-d [ --decoy ]
```
Use a shuffled decoy protein database to estimate false discovery rates. When -d is chosen, a shuffled decoy database is automatically generated and appended to the target database before database search, and false discovery rates are estimated using the target-decoy approach.
```
-e [ --error-tolerance ] <a positive integer>
```
Set the error tolerance for precursor and fragment masses in part-per-million (ppm). Default value: 15. Valid values are 5, 10, and 15 ppm. To use other values, the generating function approach (-g) should be chosen.
```
-m [ --max-shift ] <a positive number>
```
Set the maximum absolute value of the mass shift (in Dalton) of an unexpected modification. Default value: 500 Da.
```
-p [ --num-shift ] <0|1|2>
```
Set the maximum number of unexpected mass shifts in a proteoform spectrum-match. Default value: 1.
```
-t [ --spectrum-cutoff-type ] <EVALUE|FDR>
```
Set the spectrum level cutoff type for filtering proteoform spectrum-matches. Default value: EVALUE.
```
-v [ --spectrum-cutoff-value ] <a positive number>
```
Set the spectrum level cutoff value for filtering proteoform spectrum-matches. Default value: 0.01.
```
-T [ --proteoform-cutoff-type ] <EVALUE|FDR>
```
Set the proteoform level cutoff type for filtering proteoform spectrum-matches. Default value: EVALUE.
```
-V [ --proteoform-cutoff-value ] <a positive number>
```
Set the proteoform level cutoff value for filtering proteoform spectrum-matches. Default value: 0.01.
```
-g [ --generating-function ]
```
Use the generating function approach to compute p-values and E-values. It is more accurate than the default lookup table based approach, but it is much slower.
```
-r [ --num-combined-spectra ] <a positive integer>
```
Set the number of combined spectra. The parameter is set as 2 (or 3) for combining spectral pairs (or triplets) generated by the alternating fragmentation mode. Default value: 1.
```
-i [ --mod-file-name ] <a common modification file>
```
Specify a text file containing a list of common PTMs for proteoform characterization. The PTMs are used to identify and localize PTMs in reported proteoform spectrum-matches with unknown mass shifts. See example file.
```
-s [ --miscore-threshold ] <a number between 0 and 1>
```
Set the score threshold (modification identification score) for filtering results of PTM characterization. Default value: 0.45.
```
-u [ --thread-number ] <a positive number>
```
Set the number of threads used in the computation. Default value: 1.
```
-x [ --use-topfd-feature ] <a TopFD feature file with its path>
```
Specify a TopFD feature file for obtaining proteoform identifications by grouping identified proteoform spectrum-matches.
```
-l [ --skip-list ] <a text file with its path>
```
Specify a text file containing a list of spectrum scans that have been identified in previous analysis and are not needed to be searched against the protein database.
