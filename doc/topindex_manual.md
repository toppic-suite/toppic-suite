## TopIndex 
TopIndex generates index files for protein sequence databases. The index files are used in TopPIC and TopMG to speed up proteoform identification by database search.

### Minimum requirements
A computer with at least 16 GB memory and a 64-bit Linux or Windows operating system is required. 

###  Input
The input is a protein sequence database file in the FASTA format.

### Output
The output is a folder containing protein sequence index files. For example, when the input file name is proteins.fasta, the output folder is proteins.fasta_idx.

### Command line usage
To run TopIndex, open a terminal window and run the following command.
```
topindex [options] database-file-name
```

Options
```
-h [ --help ]
```
Print the help message.
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
Use a shuffled decoy protein database to estimate spectrum and proteoform-level FDRs. When -d is chosen, a shuffled decoy database is automatically generated and appended to the target database. Index files for the concatenated database are generated.
```
-e [ --mass-error-tolerance ] <a positive integer>
```
Set the error tolerance for precursor and fragment masses in ppm. Default value: 10.
```
-u [ --thread-number ] <a positive integer>
```
Set the number of threads used in the computation. Default value: 1. About 0.5 GB memory is required for each CPU thread.

### Examples

Generate index files for a protein database file proteins.fasta using default parameters.
```
topindex proteins.fasta
```

Generate index files for a protein database file proteins.fasta using carbamidomethylation as the fixed modification, and N-terminal methionine excision and N-terminal methionine acetylation as the N-terminal forms.
```
topindex -f C57 -n NME,M_ACETYLATION proteins.fasta
```

Generate index files for a protein database file proteins.fasta using a target-decoy concatenated database, a mass error tolerance of 5 ppm, and 4 CPU threads.
```
topindex -d -e 5 -u 4 proteins.fasta
```