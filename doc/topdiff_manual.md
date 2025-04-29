## TopDiff
TopDiff compares the abundances of proteoforms and finds differentially expressed proteoforms by using identifications of top-down mass spectrometry data of several protein samples.

### Minimum requirements
A computer with at least 16 GB memory and a 64-bit Linux or Windows operating system is required. 

### Input

* Proteoform identification files in the XML format, e.g., spectra_ms2_toppic_proteoform.xml

### Output

TopDiff outputs a TSV file containing proteoform identifications and their abundances in the input mass spectrum data. The default output file name is sample_diff.tsv.

### Command line usage

To run TopDiff, open a terminal window and run the following command.
```
topdiff [options] spectrum-file-names
```

Options
```
-h [ --help ]
```
Print the help message.
```
-e [ --error-tolerance ] <a positive number>
```
Set the error tolerance for mapping identified proteoforms across multiple samples (in Dalton). Default value: 1.2 Dalton.
```
-t [ --tool-name ] <toppic|topmg>
```
Specify the name of the database search tool: toppic or topmg. Default: toppic.
```
-o [ --output ] <a file name>
```
Specify the output file name. Default: sample_diff.tsv.

### Examples

Compare proteoform abundances using TopPIC identifications of two spectrum files spectra1_ms2.msalign and spectra2_ms2.msalign.
```
topdiff spectra1_ms2.msalign spectra2_ms2.msalign
```
Compare proteoform abundances using TopMG identifications of two spectrum files spectra1_ms2.msalign and spectra2_ms2.msalign.
```
topdiff -t topmg spectra1_ms2.msalign spectra2_ms2.msalign 
```