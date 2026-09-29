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

## Manuals

The command-line tools are documented in `docs/manual/`:

- [TopFD](docs/manual/topfd_manual.md) — spectral deconvolution of mzML/mzXML files and of a single peak list
- [TopDIA](docs/manual/topdia_manual.md) — pseudo-MS/MS spectra from top-down DIA data
- [TopIndex](docs/manual/topindex_manual.md) — index files for a protein database
- [TopPIC](docs/manual/toppic_manual.md) — proteoform identification and characterization
- [TopMG](docs/manual/topmg_manual.md) — identification of highly modified proteoforms with mass graphs
- [TopDiff](docs/manual/topdiff_manual.md) — differential abundance of proteoforms across samples
- [toplib_build](docs/manual/toplib_build_manual.md) — a spectral library of identified MS/MS spectra (the C++ port of TopLib's library building)
- [Proteoform annotation](docs/manual/proteoform.md) — how to read the annotated proteoform strings in TopPIC and TopMG results

The GUI tools (`topfd_gui`, `topdia_gui`, `topindex_gui`, `toppic_gui`,
`topmg_gui`, `topdiff_gui`) offer the same options as the command-line tools
and run them.

Claude Code loads the skills in `.claude/skills/` automatically when it is
started in this repository: `tools` for running the command-line tools and
reading their outputs, `compile` for building and installing the suite, and
`repo` for publishing the private repository to the public one.

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

## Citation

If you use TopPIC Suite in your research, please cite the paper of each tool
you used.

**TopFD**\
Abdul Rehman Basharat, Yong Zang, Liangliang Sun, and Xiaowen Liu. TopFD: a
proteoform feature detection tool for top-down proteomics. *Analytical
Chemistry*, 95:8189-8196, 2023.
[PubMed](https://pubmed.ncbi.nlm.nih.gov/37196155/)

**TopPIC**\
Qiang Kou, Likun Xun, and Xiaowen Liu. TopPIC: a software tool for top-down
mass spectrometry-based proteoform identification and characterization.
*Bioinformatics*, 32:3495-3497, 2016.
[PubMed](https://pubmed.ncbi.nlm.nih.gov/27423895/)

**TopMG**\
Qiang Kou, Si Wu, Nikola Tolić, Ljiljana Paša-Tolić, Yunlong Liu, and Xiaowen
Liu. A mass graph-based approach for the identification of modified
proteoforms using top-down tandem mass spectra. *Bioinformatics*,
33:1309-1316, 2017.
[PubMed](https://pubmed.ncbi.nlm.nih.gov/28453668/)

**TopDIA**\
Abdul Rehman Basharat, Xingzhao Xiong, Tian Xu, Yong Zang, Liangliang Sun, and
Xiaowen Liu. TopDIA: a software tool for top-down data-independent acquisition
proteomics. *Journal of Proteome Research*, 24:55-64, 2025.
[PubMed](https://pubmed.ncbi.nlm.nih.gov/39641251/)

**FLASHDeconv**\
TopFD deconvolutes MS/MS spectra with the FLASHDeconv algorithm when
`--ms2-flashdeconv` is given. If you use that option, please also cite:\
Kyowon Jeong, Jihyung Kim, Manasi Gaikwad, Siti Nurul Hidayah, Laura Heikaus,
Hartmut Schlüter, and Oliver Kohlbacher. FLASHDeconv: ultrafast, high-quality
feature deconvolution for top-down proteomics. *Cell Systems*, 10:213-218,
2020.
[PubMed](https://pubmed.ncbi.nlm.nih.gov/32078799/)
