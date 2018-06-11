To run TopMG, a computer with at least 4 GB memory and a 64-bit Linux or Windows operating system is required. TopPIC provides a command line interface and a graphical user interface (GUI) for both Linux and Windows users.

# TopMG

## Input

* A protein database file in the FASTA format
* A mass spectrum data file in the msalign format
* A text file to specify variable PTMs for proteoform identification
* A text file to specify fixed PTMs (optional)
* A text file containing MS1 feature information (optional)

## Output

TopMG outputs two tab delimited text files and a collection of html files for identified proteoforms. For example, when the input mass spectrum data file is spectra_ms2.msalign, the output includes:

* __spectra_ms2.OUTPUT_TABLE__: a tab delimited text file containing a list of identified proteoform spectrum-matches at the spectrum level.
* __spectra_ms2.FORM_OUTPUT_TABLE__: a tab delimited text file containing a list of identified proteoform spectrum-matches at the proteoform level.
* __spectra_ms2_prsm_cutoff_html__: a folder containing html files for the annotation of identified proteoform spectrum-matches using the spectrum level E-value or FDR cutoff.
* __spectra_ms2_proteoform_cutoff_html__: a folder containing html files for the annotation of identified proteoform spectrum-matches using the proteoform level E-value or FDR cutoff.

## Command line usage

To run TopMG, open a console and run the following command.
```
topmg [options] database-file-name spectrum-file-name
```
Options
```
-h [ --help ]
```
Print the help message.

```
-a [ --activation ] <CID|HCD|ETD|UVPD|FILE>
```
Fragmentation method of tandem mass spectra. When FILE is used, fragmentation methods of spectra are given in the input spectral data file. Default value: FILE.

```
-f [ --fixed-mod ] <C57|C58|a fixed modification file>
```
Fixed modifications. Three available options: C57, C58, or the name of a text file containing the information of fixed modifications. When C57 is selected, carbamidomethylation on cysteine is the only fixed modification. When C58 is selected, carboxymethylation on cysteine is the only fixed modification.

```
-n [ --n-terminal-form ] <a list of allowed N-terminal forms>
```
N-terminal forms of proteins. Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, and M_ACETYLATION. NONE stands for no modifications, NME for N-terminal methionine excision, NME_ACETYLATION for N-terminal acetylation after the initiator methionine is removed, and M_ACETYLATION for N-terminal methionine acetylation. When multiple forms are allowed, they are separated by commas. Default value: NONE,NME,NME_ACETYLATION,M_ACETYLATION.

```
-d [ --decoy ]
```
Use a decoy protein database to estimate false discovery rates.

```
-e [ --error-tolerance ] <a positive integer>
```
Error tolerance for precursor and fragment masses in PPM. Default value: 15.

```
-m [ --max-shift ] <a positive number>
```
Maximum absolute value of the mass shift (in Dalton). Default value: 500.

```
-t [ --spectrum-cutoff-type ] <EVALUE|FDR>
```
Spectrum-level cutoff type for filtering identified proteoform spectrum-matches. Default value: EVALUE.

```
-v [ --spectrum-cutoff-value ] <a positive number>
```
Spectrum-level cutoff value for filtering identified proteoform spectrum-matches. Default value: 0.01.

```
-T [ --proteoform-cutoff-type ] <EVALUE|FDR>
```
Proteoform-level cutoff type for filtering identified proteoform spectrum-matches. Default value: EVALUE.

```
-V [ --proteoform-cutoff-value ] <a positive number>
```
Proteoform-level cutoff value for filtering identified proteoform spectrum-matches. Default value: 0.01.

```
-i [ --mod-file-name ] <a common modification file>
```
Specify a text file containing the information of common PTMs for constructing proteoform graphs.

```
-u [ --thread-number ] <a positive integer>
```
Number of threads used in the computation. Default value: 1.

```
-x [ --no-topfd-feature ]
```
No TopFD feature file for proteoform identification.

```
-l [ --skip-list ] <a text file with its path>
```
The scans in this file will be skipped.

```
-o [ --output ] <a filename with its path>
```
The output file name for the combined results. Default: combined.

```
-j [ --proteo-graph-dis ] <a positive number>
```
Gap in constructing proteoform graph. Default value: 40.

```
-G [ --var-ptm-in-gap ] <a positive number>
```
Maximum number of variable PTMs in a proteform graph gap. Default value: 5.

```
-D [ --use-asf-diagonal ]
```
Use the ASF-DIAGONAL method for protein filtering.

```
-P [ --var-ptm ] <a positive number>
```
Maximum number of variable PTMs. Default value: 5.
