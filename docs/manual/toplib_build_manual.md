# toplib_build manual

toplib_build builds a top-down spectral library: it clusters identified
MS/MS spectra, makes one representative spectrum per cluster, gives it the
cluster's best identification and a decoy spectrum, and writes the targets
and decoys into a SQLite database. It is a C++ port of the library building
of TopLib (`src/ms_library_building.py`) and writes the same table, with two
differences: spectra identified as different proteoforms are never put in
one cluster (step 2 in section 2), and the PrSM table's columns are also
found under the names TopPIC 1.9 and 2.0 give them (section 1).

```sh
toplib_build [options] <msalign file> [<PrSM file>]
```

`toplib_build -h` prints the option list.

## 1. Input

Either of two inputs:

1. **A TopFD msalign file and a TopPIC PrSM table** (two arguments), the
   input of TopLib's script, for example `sample_ms2.msalign` and
   `sample_ms2_toppic_prsm_single.tsv`. From each spectrum of the msalign
   file toplib_build reads `FILE_NAME`, `SPECTRUM_ID`, `SCANS`,
   `RETENTION_TIME`, `PRECURSOR_MASS`, `PRECURSOR_INTENSITY`,
   `PRECURSOR_CHARGE`, `PRECURSOR_FEATURE_ID` (the first value where a
   spectrum lists several precursors) and, if present, `TITLE`, `LEVEL`,
   `MS_ONE_ID`, `MS_ONE_SCAN`, `PRECURSOR_WINDOW_BEGIN`,
   `PRECURSOR_WINDOW_END`, `ACTIVATION` and `PRECURSOR_MZ`; each field is the
   first header line that contains its name. The PrSM table is joined to the
   spectra by `Spectrum ID`. Its header is its first line with a
   `Spectrum ID` column, so TopPIC's parameter block above the header is
   skipped. It must have the columns `Proteoform ID`, `Protein accession`,
   `Protein description`, `First residue`, `Last residue`,
   `Database protein sequence`, `Proteoform`, `E-value`,
   `Spectrum-level Q-value` and `Proteoform-level Q-value`; these columns
   are copied when present: `#unexpected modifications`,
   `unexpected modifications`, `Fixed PTMs`, `Adjusted precursor mass`,
   `Prsm ID`, `Fragmentation`, `#peaks`, `Proteoform intensity`,
   `Feature intensity`, `Feature score`, `Feature apex time`,
   `#Protein hits`, `Proteoform mass`, `Protein N-terminal form`,
   `#variable PTMs`, `variable PTMs`, `#matched peaks`,
   `#matched fragment ions`, `Special amino acids` and `MIScore`. These are
   TopLib's names, those of TopPIC 1.7's table; TopPIC 1.9 and 2.0 renamed
   eight of them, and those names are accepted as well:

   | TopLib (TopPIC 1.7) | TopPIC 1.9 | TopPIC 2.0 |
   |---|---|---|
   | `#unexpected modifications` | `# Unexpected modifications` | `# Unexpected modifications` |
   | `unexpected modifications` | `Unexpected modifications` | `Unexpected modifications` |
   | `#peaks` | `# Masses` | `# Masses` |
   | `#Protein hits` | `# Protein hits` | `# Protein hits` |
   | `#variable PTMs` | `# Variable PTMs` | `# Variable PTMs` |
   | `variable PTMs` | `Variable PTMs` | `Variable PTMs` |
   | `#matched peaks` | `# Matched masses` | `# Matched exp fragments` |
   | `#matched fragment ions` | `# Matched fragments` | `# Matched theo fragments` |

   Spectra without a row in the table are unidentified.
2. **A library msalign file** (one argument): an msalign file whose spectra
   carry their identification in the header (`PROTEIN_ACCESSION`,
   `DATABASE_SEQUENCE`, `FIRST_RESIDUE_POSITION`, `PROTEOFORM`,
   `FIXED_MODIFICATIONS`, `UNEXPECTED_MODIFICATIONS`, `E_VALUE`) and name
   the spectrum fields `MS2_ID`, `MS2_SCAN`, `MS2_RETENTION_TIME`, `MS1_ID`,
   `MS1_SCAN`, `PRECURSOR_MONOISOTOPIC_MZ`, `PRECURSOR_MONOISOTOPIC_MASS`
   and `MZML_FILE_NAME`. The fields map onto the columns of the first input;
   the proteoform id is the index of the spectrum's (accession, proteoform)
   pair in the file, the last residue is the first residue plus the number
   of residues of the proteoform minus one, and the protein description and
   the q-values are empty. The spectrum id is `MS2_ID`.

In both, spectra without masses or with a precursor mass of 0 are skipped,
and masses above the precursor mass are left out.

## 2. What toplib_build does

1. **Reads the spectra** (`Reading spectra`) and, with `-c`, keeps the
   spectra whose first precursor has at least that fraction of the summed
   intensity of the spectrum's precursors (`PPIR filter`).
2. **Clusters the spectra** (`Clustering spectra and building
   representatives`). Spectra with fewer than two masses or a precursor
   intensity of 0 are left out. The others are grouped by precursor charge,
   each group is clustered by precursor mass (complete linkage, clusters
   within 2.2 Da), and each of those clusters by spectral similarity
   (average linkage, cosine distance at most 0.7). The cosine distance uses
   the 50 most intense masses of each spectrum with log2 intensities
   scaled to unit length; two masses of equal charge match when they differ
   by at most 10 ppm, or by 1.00235 Da within 10 ppm. Two spectra
   identified as different proteoforms (a different proteoform id or
   protein accession) are never put in one cluster, however similar their
   spectra; unidentified spectra can join any cluster. (TopPIC 1.9 could
   give PrSMs of different proteins one proteoform id; 2.0 does not.) (TopLib's script lets them share a cluster, whose
   representative then keeps only one of their identifications, so the
   other proteoforms can be left without a representative.)
3. **Makes one representative per cluster.** `average` (the default)
   merges the members' spectra (all masses): masses of equal charge within
   10 ppm are combined into their intensity-weighted mean mass with the
   summed intensity. Its precursor mass is the intensity-weighted mean of
   the members whose integer precursor mass differs from the first
   member's by the most common difference, and its precursor intensity their
   mean intensity. `single` takes the member with the lowest E-value (the
   member with the most masses when no member is identified). The
   representative keeps its 50 most intense masses, sorted by mass.
4. **Identifies each representative** with the member of lowest E-value.
   Of the clusters that share a proteoform (id and accession) and a
   precursor charge, only the one with the lowest E-value keeps its
   identification; representatives
   without an identification are not written.
5. **Makes a decoy** of each representative: every mass is shifted by the
   mass difference of two randomly chosen amino acid residues (redrawn
   while the shifted mass is negative, then while it is above the precursor
   mass), from
   Python's random number generator seeded with the representative's index,
   so the decoys are the same as TopLib's.
6. **Writes the library** (`Writing the library`).

The console prints the number of spectra read and identified, the number of
precursor clusters and clusters, how many clusters have no identification or
lose it, and the number of targets and decoys written.

## 3. Output

A SQLite database (default `TopLib/lib_spectra_ms2.db`, the directory is
created) with the table `target_decoy_spectra_representatives`; an existing
table of that name is replaced and other tables are kept. Each row is a
target (`flag` 1) or a decoy (`flag` 0) representative; the targets come
first, and a decoy has its target's columns apart from the masses.

| Column | Meaning |
|---|---|
| `representative_id` | The cluster's number (from 1), shared by the target and the decoy; clusters without an identification leave gaps. |
| `precursor_mass`, `precursor_intensity`, `precursor_charge` | The representative's precursor. |
| `mass_0` ... `mass_49` | The masses by increasing mass (the decoy's shifted masses in the same positions); empty after the last mass. |
| `intensity_0` ... `intensity_49`, `charge_0` ... `charge_49` | Their intensities and charges. |
| `norm_intensity_0` ... `norm_intensity_49` | log2 of the intensities, scaled to unit length. |
| `flag` | 1 target, 0 decoy. |
| `file_id` | 1. |
| `file_name`, `spectrum_id`, `title`, `scan`, `retention_time`, `ms_level`, `ms_one_id`, `ms_one_scan`, `precursor_window_begin`, `precursor_window_end`, `activation`, `precursor_mz`, `precursor_feature_id` | The spectrum of the identification (the cluster member with the lowest E-value). The optional fields are columns only when some spectrum has them. |
| `proteoform_id`, `protein_accession`, `protein_description`, `first_residue`, `last_residue`, `protein_sequence`, `proteoform`, `e_value`, `spectrum_level_q_value`, `proteoform_level_q_value`, then the copied optional columns (`unexpected_modifications`, `fixed_ptms`, ...) | Its identification. |
| `cluster_id` | The cluster's number from 0. |

The column types are those TopLib's pandas tables get: `charge_k` is
INTEGER when every representative has more than k masses and REAL
otherwise, the identification columns are REAL when all their values are
numbers and TEXT otherwise, and `cluster_id` is REAL.

## 4. Options

| Option | Default | Meaning |
|---|---|---|
| `-r`, `--representative <average\|single>` | average | The representative spectrum of a cluster: the merged spectrum of its members, or its member with the lowest E-value. |
| `-c`, `--ppir-cutoff <0..1>` | 0 | Keep only the spectra whose first precursor has at least this fraction of the summed intensity of the spectrum's precursors (the precursor-to-precursor intensity ratio); 0 keeps all spectra. |
| `-o`, `--output <file>` | TopLib/lib_spectra_ms2.db | The library database. |

## 5. Examples

A library from a TopFD msalign file and TopPIC's PrSM table, the
representatives merged from the cluster members:

```sh
toplib_build sample_ms2.msalign sample_ms2_toppic_prsm_single.tsv
```

writes `TopLib/lib_spectra_ms2.db` in the current directory. The same with
the best-scoring spectrum of each cluster, from spectra with a PPIR of at
least 0.5:

```sh
toplib_build -r single -c 0.5 -o single_lib.db sample_ms2.msalign sample_ms2_toppic_prsm_single.tsv
```

A library from a library msalign file:

```sh
toplib_build -o lib.db sw480_rep1_lib_v1.2.1.msalign
```
