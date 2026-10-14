# toplib_search manual

toplib_search identifies top-down MS/MS spectra by searching them against a
spectral library built by toplib_build (or TopLib). Each query spectrum is
compared with the library spectra of a matching precursor mass, it takes
the most similar one, and the target matches that pass a false discovery
rate (FDR) cutoff are written to a table. It is a C++ port of the library
search of TopLib (`src/ms_library_query.py`) and writes the same table.

```sh
toplib_search [options] <library file> <msalign file>
```

`toplib_search -h` prints the option list.

## 1. Input

1. **A library**: a SQLite database with the table
   `target_decoy_spectra_representatives`, as toplib_build writes it
   (default name `TopLib/lib_spectra_ms2.db`). Its rows with a proteoform
   are searched, targets (`flag` 1) and decoys (`flag` 0). The file is
   used as given (TopLib's script looks for it under `./TopLib/`).
2. **An msalign file** of deconvoluted MS/MS spectra, as TopFD writes it,
   for example `sample_ms2.msalign`. From each spectrum toplib_search reads
   `FILE_NAME`, `SPECTRUM_ID`, `PRECURSOR_MASS`, `PRECURSOR_INTENSITY` and
   `PRECURSOR_CHARGE` (the first value where a spectrum lists several
   precursors; each field is the first header line that contains its name)
   and the masses. Spectra without masses or with a precursor mass of 0 are
   skipped, masses above the precursor mass are left out, and spectra with
   a precursor intensity of 0 or fewer than two masses are not searched.

## 2. What toplib_search does

1. **Reads the library** (`Reading the library`) and **the spectra**
   (`Reading spectra`).
2. **Searches the library** (`Searching the library`). For each query
   spectrum:
   - The candidates are the library spectra whose precursor mass matches
     the query's: by default within 10 ppm of the library spectrum's
     precursor mass, or one isotope shift (1.00228 Da) away within 10 ppm
     (`-E`); with `-T Da`, within 2.2 Da. With `-c yes` the precursor
     charge must be the same too.
   - The query's 50 most intense masses are sorted by mass and their log2
     intensities scaled to unit length, as the library spectra's are.
   - The cosine distance to each candidate is one minus the sum of the
     intensity products of the mass pairs of equal charge that differ by
     at most 10 ppm (`-e`, in ppm of the library mass) or by one isotope
     shift within it; when the sum passes 1 it stops, without the last
     product.
   - The query takes the first candidate of the smallest distance. The
     console reports how many queries got a library spectrum and how many
     of those are targets and decoys.
3. **FDR** (`FDR - cosine distance cutoff`). For each distance cutoff 0.1,
   0.101, ..., 0.999 the FDR is the decoy fraction of the matches at or
   below it; the cutoff whose FDR is closest to 1% (above or below) is
   chosen. The console prints it with its target and decoy numbers.
4. **Writes the identifications** (`Writing the identifications`): the
   target matches at the cutoff with a similarity (one minus the distance)
   of at least 0.3.

## 3. Output

A tab-separated table (default `toplib_output/query_res.tsv`, the directory
is created), one row per identified spectrum, by spectrum id:

| Column | Meaning |
|---|---|
| `spectrum_id`, `file_name` | The query spectrum. |
| `representative_id` | The library spectrum. |
| `similarity` | One minus the cosine distance. |
| `file_name(lib)`, `spectrum_id(lib)`, `scan(lib)`, `retention_time(lib)` | The spectrum of the library spectrum's identification. |
| `proteoform_id(lib)`, `protein_accession(lib)`, `precursor_mass(lib)`, `precursor_charge(lib)`, `precursor_feature_id(lib)`, `protein_description(lib)`, `protein_sequence(lib)`, `first_residue(lib)`, `last_residue(lib)`, `proteoform(lib)`, `e_value(lib)`, `spectrum_level_q_value(lib)`, `proteoform_level_q_value(lib)` | The library spectrum's precursor and identification. |

The values are written as TopLib writes them: the numbers of a library
column that holds real numbers or has an empty value are written as real
numbers (`924.0`), and empty values as empty cells.

## 4. Options

| Option | Default | Meaning |
|---|---|---|
| `-T`, `--precursor-error-type <ppm\|Da>` | ppm | The unit of the precursor mass error tolerance. |
| `-E`, `--precursor-error <number>` | 10 ppm or 2.2 Da | The precursor mass error tolerance. In ppm (of the library spectrum's precursor mass), a library spectrum matches when its precursor mass is within the tolerance of the query's or one isotope shift (1.00228 Da) away from it within the tolerance; in Da, when it is within the tolerance. |
| `-e`, `--fragment-error <number>` | 10 | The fragment mass error tolerance in ppm. |
| `-c`, `--charge-state <yes\|no>` | no | Whether a library spectrum must have the query's precursor charge. |
| `-o`, `--output <file>` | toplib_output/query_res.tsv | The table of identifications. |

## 5. Examples

Search TopFD's spectra against a library built by toplib_build, with the
default 10 ppm precursor and fragment tolerances:

```sh
toplib_build sample_ms2.msalign sample_ms2_toppic_prsm_single.tsv
toplib_search TopLib/lib_spectra_ms2.db query_ms2.msalign
```

writes `toplib_output/query_res.tsv`. The same with a 2.2 Da precursor
window, the same precursor charge required and another output file:

```sh
toplib_search -T Da -c yes -o query_da.tsv TopLib/lib_spectra_ms2.db query_ms2.msalign
```
