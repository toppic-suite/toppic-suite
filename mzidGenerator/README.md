# mzidGenerator

Compatible with python2 and python3

## Installation
```sh
#only if Python is not installed already
sudo apt-get install python3.6
            
#use pip2 if using python2
pip3 install lxml
pip3 install six
pip3 install sqlalchemy
pip3 install numpy
pip3 install pyteomics
```
## Usage
```sh
python write_mzIdent.py <path-to-prsm-tsv-file> <path-to-fasta-or-fasta_target_decoy-file> <optional: path-to-fixed-ptm-file> <optional: path-to-common-PTM-file> <optional: path-to-variable-PTM-file>
```

### Examples
```sh
#use None if optional parameters are not used
#example: using tsv file from TopPIC
python write_mzIdent.py st_1_ms2_toppic_prsm.tsv uniprot-st.fasta fixed_PTM.txt common_PTM.txt None

#example: using tsv file from TopMG
python write_mzIdent.py st_1_ms2_toppic_prsm.tsv uniprot-st.fasta fixed_PTM.txt None variable_PTM.txt
```
