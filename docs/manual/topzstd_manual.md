# topzstd manual

`topzstd` compresses the SQLite databases of TopFD, TopDIA, TopPIC and
TopMG, or decompresses them. topfd, topdia, toppic and topmg compress their
databases by default; `topzstd` compresses the databases of earlier
versions, and gives the rows back (`-d`) for programs that read SQLite
without the sqlite-zstd extension.

```
topzstd [options] database-file-name ...
```

## 1. The compressed format

The tables of peaks and envelopes (`ms1_peak`, `ms2_peak`, `ms1_env`,
`ms1_env_peak`, `ms2_env`, `ms2_env_peak`, whichever the database has) are
stored as one row per spectrum, with one array per column, and the arrays
are compressed with zstd:

- An m/z or intensity column is stored as 32-bit floats when every one of
  its values is exactly a 32-bit float (the peaks of TopFD), otherwise as
  64-bit floats; integers take the smallest integer type that holds them.
  Ids that number the rows of each spectrum (or of each envelope) 0, 1, ...
  are not stored. Nothing is lost: decompression gives every value back.
- Each column is compressed with a zstd dictionary of its own, trained on
  all its arrays; a column too small for a dictionary (under about 25 kB) is
  stored uncompressed.
- Every other table (the spectra, the 3D levels' `CONFIG`, TopPIC's and
  TopMG's identification tables, the search database) is copied as it is.

The format is that of the [sqlite-zstd](https://github.com/phiresky/sqlite-zstd)
extension (each compressed table is a view `<table>_blob` that the
extension decompresses), so TopMSV reads the files with the extension
loaded. toppic and topmg read a compressed TopFD database themselves.
A compressed database of a TopFD run is about a third of its size (for
example 843 MB -> 295 MB); a TopPIC database, which also holds the search
database, a half to four fifths.

## 2. Input

One or more SQLite databases of TopFD (`sample.sqlite`), TopDIA, TopPIC
(`sample_toppic.sqlite`) or TopMG (`sample_topmg.sqlite`). Every input is
checked before any is changed: an input that is not an SQLite database, a
compressed input (or with `-d` an uncompressed one) stops the run with an
error.

## 3. Output

Each input is replaced by its compressed (or with `-d` decompressed) form,
unless `-o` names another file. The new database is written under the
output name with `.tmp` appended and renamed when complete, so an
interrupted run leaves the input as it was. For each file, `topzstd`
prints the sizes before and after and the time.

A decompressed database has every table and value of the database before
compression (the same CREATE statements); the rows of each spectrum are in
their order, the spectra in the order of their ids.

## 4. Options

| Option | Default | Description |
|---|---|---|
| `-h`, `--help` | | Print the help message. |
| `-d`, `--decompress` | off | Decompress: write the database with its peak and envelope tables back as rows, as it was before compression. |
| `-l`, `--level <1-22>` | 3 | The zstd compression level. |
| `-o`, `--output <a file name>` | | Write the result to this file instead of replacing the input file (one input file only). |

## 5. Examples

```
# compress the databases of an earlier TopFD and TopPIC run in place
topzstd sample.sqlite sample_toppic.sqlite

# a decompressed copy, for a program without the sqlite-zstd extension
topzstd -d sample.sqlite -o sample_rows.sqlite
```
