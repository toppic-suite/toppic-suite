To run TopFD, a computer with at least 4 GB memory and a 64-bit Linux or Windows operating system is required. TopFD provides a command line interface and a graphical user interface (GUI) for both Linux and Windows users.

# TopFD

## Input

The input of TopFD is mzML or mzXML top-down mass spectrometry data files. Raw mass spectral data generated from various mass spectrometers can be converted to mzML or mzXML files using msconvert.

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
topfd [options] spectrum-file-names
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
