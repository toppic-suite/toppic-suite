## TopFD
TopFD (Top-down mass spectral Feature Detection) is a software tool for top-down spectral deconvolution, which groups top-down mass spectral peaks into isotopic envelopes and converts isotopic envelopes to monoisotopic neutral masses. In addition, it extracts proteoform features from MS1 spectra. 

### Minimum requirements
A computer with at least 16 GB memory and a 64-bit Linux or Windows operating system is required. 

### Input
The input of TopFD is mzML or mzXML top-down mass spectrometry data files. Raw mass spectral data generated from various mass spectrometers can be converted to mzML or mzXML files using [MSConvertGUI](https://proteowizard.sourceforge.io/tools/tools_base.html).


### Output
TopFD outputs two LC-MS feature text files with a file extension `feature`, one LC-MS feature file with a file extension `xml`, and two deconvoluted mass spectral data files in the msalign format with a file extension `msalign`, which is similar to the MGF file format. In addition, TopFD creates a folder containing JavaScript files for spectral visualization.

For example, when the input file name is spectra.mzML, the output includes: 

* **spectra_ms1.feature**: a feature file containing LC-MS features.
* **spectra_ms2.feature**: a feature file containing MS/MS scan IDs and their corresponding LC-MS feature IDs.
* **spectra_feature.xml**: a feature file containing LC-MS features in the xml format.
* **spectra_ms1.msalign**: a list of deconvoluted MS1 spectra.
* **spectra_ms2.msalign**: a list of deconvoluted MS/MS spectra.
* **spectra_html**: a folder containing JavaScript files for MS1 and MS/MS spectral visualization.


### Command line usage
To run TopFD, open a console and run the following command.
```
topfd [options] spectrum-file-names
```

Options

```
-h [ --help ]
```
Print the help message.

```
-a [ --activation ] <CID|ETD|HCD|MPD|UVPD|FILE>
```
Set the fragmentation method(s) of MS/MS spectra. When "FILE" is selected, the fragmentation methods of spectra are given in the input spectrum data file. Default value: FILE.

```
-c [ --max-charge ] <a positive integer>
```
Set the maximum charge state of precursor and fragment ions. The default value is 30.

```
-m [ --max-mass ] <a positive number>
```

Set the maximum monoisotopic mass of precursor and fragment ions. The default value is 50,000 Dalton.
```
-e [ --mz-error ] <a positive number>
```
Set the error tolerance of m/z values of spectral peaks. The default value is 0.02 m/z.
```
-r [ --ms-one-sn-ratio ] <a positive number>
```
Set the signal/noise ratio for MS1 spectra. The default value is 3.
```
-s [ --ms-two-sn-ratio ] <a positive number>
```
Set the signal/noise ratio for MS/MS spectra. The default value is 1.
```
-o [ --missing-level-one ]
```
Specify that the input file does not contain MS1 spectra.
```
-n [ --msdeconv ]
```
Use the MS-Deconv score (see [paper](https://pubmed.ncbi.nlm.nih.gov/20855543/)) to rank isotopic envelopes. If -n is not selected, the default EnvCNN score (see [paper](https://pubmed.ncbi.nlm.nih.gov/32356965/)) is used to rank isotopic envelopes.
```
-w [ --precursor-window ] <a positive number>
```
Set the precursor isolation window size. The default value is 3.0 m/z. When the input file contains the information of precursor windows, the parameter will be ignored.
```
-t [ --ecscore-cutoff ] <a number in [0, 1]>
```
Set the ECScore cutoff value for proteoform features. Default value is 0.5.
```
-b [ --min-scan-number ] <1|2|3>
```
The minimum number of MS1 scans in which a proteoform feature is detected. The default value is 3.
```
-i [ --single-scan-noise ]
```
Use the noise intensity levels in single MS1 scans to filter out low intensity peaks in proteoform feature detection. The default method is to use the noise intensity level of the whole LC-MS map to filter out low intensity peaks.
```
-f [ --additional-feature-search ]
```
Perform additional proteoform feature search in the LC-MS map for MS/MS scans that do not have detected proteoform features in their precursor isolation windows. In the additional search, the signal noise ratio is set to 0, the mininum scan number is set to 1, and the ecscore cutoff is set to 0.
```
-d [ --disable-frag-num-filtering ]
```
Skip filtering fragment ion envelopes in an MS/MS spectrum based on the estimated number of fragment ions. For CID or HCD MS/MS spectra, the expected numbers of b- and y-ions are estimated and used for filtering. For ETD spectra, the estimated numbers of c- and z•-ions guide the filtering process

```
-u [ --thread-number ] <a positive integer>
```
Number of CPU threads used in spectral deconvolution. Default value: 1.
```
-g [ --skip-html-folder ]
```
Skip the generation of HTML files for visualization.


### Examples
Deconvolute a centroid data file spectra.mzML and output five files: spectra_ms1.feature, spectra_ms2.feature, spectra_feature.xml, spectra_ms1.msalign, and spectra_ms2.msalign.
```
topfd spectra.mzML
```
Deconvolute two centroid data files spectra1.mzML and spectra2.mzML and output five files for each input data file.
```
topfd spectra1.mzML spectra2.mzML
```
Deconvolute all centroid data files in the current folder.
```
topfd *.mzML
```
Deconvolute a centroid data file spectra.mzML and skip the final filtering and skip the generatation of the HTML folder for visualization.
```
topfd -d -g spectra.mzML
```
Deconvolute a centroid data file spectra.mzML using 4 CPU threads and MS-deconv score.
```
topfd -u 4 -n spectra.mzML
```
Deconvolute a centroid data file spectra.mzML that does not contain MS1 spectra.
```
topfd -o spectra.mzML
```
Deconvolute a centroid data file spectra.mzML. In proteoform feature identification, each proteoform feature is required to be detected in at least one MS1 scan and the ECScore cutoff is set to 0.2. This settings will increase the number of reported proteoform features.
```
topfd -t 0.2 -b 1 spectra.mzML
```
Deconvolute a centroid data file spectra.mzML with a signal/noise ratio 2 for MS1 spectra.
```
topfd -r 2 spectra.mzML
```
Deconvolute a centroid data file spectra.mzML with the following settings: the maximum charge state: 50, the maximum mass: 30,000 Dalton, and the signal/noise ratio for MS/MS spectra: 2.
```
topfd -c 50 -m 30000 -s 2 spectra.mzML 
```