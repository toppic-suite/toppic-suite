## TopPIC Suite

For tutorials, publications and binary downloads, please visit
https://www.toppic.org/software/toppic/. The user manuals are in this
repository (see [Manuals](#manuals) below).

TopPIC Suite consists of six software tools for the analysis of top-down mass spectrometry-based proteomics data. 

* **TopFD** (Top-down mass spectral Feature Detection) is a software tool for top-down mass spectral deconvolution. It groups top-down spectral peaks into isotopic envelopes and converts isotopic envelopes to monoisotopic neutral masses. In addition, it extracts proteoform features from LC-MS or CE-MS data.

* **TopDIA** is a tool for demultiplexing top-down data independent acquisition mass spectrometry (TD-DIA-MS) data. It processes TD-DIA-MS data to generate demultiplexed pseudo-MS/MS spectra, which are subsequently searched against a protein database to identify proteoforms.

* **TopIndex** (Top-down protein sequence database Indexing) generates index files for protein sequence databases. The index files are used in TopPIC and TopMG to speed up proteoform identification by database search.

* **TopPIC** (Top-down mass spectrometry-based Proteoform Identification and Characterization) identifies and characterizes proteoforms at the proteome level by searching top-down tandem mass spectra against a protein sequence database. It efficiently identifies proteoforms with post-translational modificatons (PTMs) and unexpected alterations, such as mutations, accurately estimates the statistical significance of identifications, and characterizes reported proteoforms with unknown mass shifts. It uses several techniques, such as indexes, spectral alignment, generating function methods, and the modification identification score (MIScore), to increase the speed, sensitivity, and accuracy.

* **TopMG** (Top-down mass spectrometry-based proteoform identification using Mass Graphs) is a software tool for identifying highly modified proteoforms by searching top-down tandem mass spectra against a protein sequence database. It is capable of identifying proteoforms with multiple variable PTMs and unexpected alterations, such as histone proteoforms and phosphorylated ones. It uses mass graphs, which efficiently represent candidate proteoforms with multiple variable PTMs, to increase the speed and sensitivity in proteoform identification. In addition, approximate spectrum-based filtering methods are employed for protein sequence filtering, and a Markov chain Monte Carlo method (TopMCMC) is used for estimating the statistical significance of identifications.

* **TopDiff** (Top-down mass spectrometry-based identification of Differentially expressed proteoforms) compares the abundances of proteoforms and finds differentially expressed proteoforms by using identifications of top-down mass spectrometry data of several protein samples.

* **TopDIA** is a software tool for top-down data-independent-acquistion mass spectrometry (TD-DIA-MS) data analysis. It generates demultiplexed pseudo MS/MS spectra from TD-DIA-MS data, which are then searched against a protein sequence database using TopPIC or TopMG for proteoform identification.

## Manuals

The command-line tools are documented in `docs/manual/`:

- [TopFD](docs/manual/topfd_manual.md) — spectral deconvolution of mzML/mzXML files and of a single peak list
- [TopDIA](docs/manual/topdia_manual.md) — pseudo-MS/MS spectra from top-down DIA data
- [TopIndex](docs/manual/topindex_manual.md) — index files for a protein database
- [TopPIC](docs/manual/toppic_manual.md) — proteoform identification and characterization
- [TopMG](docs/manual/topmg_manual.md) — identification of highly modified proteoforms with mass graphs
- [TopDiff](docs/manual/topdiff_manual.md) — differential abundance of proteoforms across samples
- [Proteoform annotation](docs/manual/proteoform.md) — how to read the annotated proteoform strings in TopPIC and TopMG results

The GUI tools (`topfd_gui`, `topdia_gui`, `topindex_gui`, `toppic_gui`,
`topmg_gui`, `topdiff_gui`) offer the same options as the command-line tools
and run them.

**TopMSV** (Top-down Mass Spectrometry Visualization) is the web-based
viewer for TopPIC Suite results. It reads the SQLite database written by
TopFD and TopPIC and shows the MS1 and MS/MS spectra with their
deconvoluted envelopes, the LC-MS map with the detected proteoform
features, and the identified proteoform-spectrum matches with their
annotated fragment ions. It is distributed separately; see the
[TopMSV manual](https://toppic.org/software/topmsv/manual.html).

## Building TopPIC Suite

* [Ubuntu Linux](docs/compile/ubuntu_build.md)
* [Red Hat Enterprise Linux](docs/compile/redhat_build.md)
* [Windows](docs/compile/windows_build.md)
* [macOS](docs/compile/macos_build.md)
