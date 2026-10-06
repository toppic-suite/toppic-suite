# Proteoform annotation in TopPIC and TopMG results

TopPIC and TopMG report every identified proteoform as an annotated
sequence string, for example

```text
MKVRASVK(KL(C)[Carbamidomethyl]RN(C)[Carbamidomethyl]K)[-116.0905]IVKRDGVIRVI(C)[Carbamidomethyl]SAEPKHKQRQG
```

This page explains how to read the string: what the parentheses, brackets,
nesting and leading `[...]-` mean, which modifications are shown by name
and which as a mass, and where the names and masses come from. Section 6
describes the `ProForma` column, the same proteoform in the ProForma
notation of the HUPO-PSI.

## 1. Where the annotation appears

- The `Proteoform` column of the PrSM, proteoform and protein tables
  (`*_prsm.tsv`, `*_proteoform.tsv`, `*_protein.tsv` and their `_single`
  versions). The residues flanking the proteoform in the protein are in the
  separate `Previous amino acid` and `Next amino acid` columns; each is `-`
  when the proteoform reaches that end of the protein.
- The `proteo_match_seq` element of the result XML files and the
  `proteoform` column of the identification tables in the SQLite database,
  which hold the same string.
- Some other outputs (visualization tools and older TopPIC releases) join
  the three columns into one string, `previous.proteoform.next`. A dot with
  nothing beside it then means the proteoform reaches the protein terminus:
  `.MKVR...QG.` is a complete protein.

## 2. Reading the string

The plain letters are the amino acids of the proteoform, from its first
residue (`First residue` column) to its last (`Last residue`), as in the
database sequence. Modifications are annotated on top of them:

| Notation | Meaning |
|---|---|
| `(X)[Name]` | Residue `X` carries the modification `Name`. |
| `(XYZ...)[label]` | The mass shift or modification `label` is located somewhere on the residues in parentheses; the fragment ions do not determine which one. |
| `(...)[+123.4567]`, `(...)[-116.0905]` | An unexpected mass shift, given as a mass in Da (see section 4). |
| `[Acetyl]-M...` | Acetylation of the protein N-terminus, written as a prefix before the first residue. |

Every annotation is written as a pair of parentheses around the residues it
applies to, immediately followed by its label in brackets. The label is
either a modification name or a signed mass.

**Nested parentheses.** Annotations are added in layers: first the
modifications present in the input, then fixed modifications, then protein
N-terminal modifications, then variable modifications, and finally the
unexpected mass shifts. When an unexpected shift spans residues that carry a
named modification, the named modification keeps its own `(X)[Name]` and the
shift's parentheses enclose it, giving nested parentheses. In the example
above,

```text
(KL(C)[Carbamidomethyl]RN(C)[Carbamidomethyl]K)[-116.0905]
```

means: both cysteines are carbamidomethylated (a fixed modification), and an
additional unexpected mass shift of -116.0905 Da is located somewhere on the
seven residues `KLCRNCK`. The inner annotations are per residue; the outer
one is a range.

**Position ranges.** The residues inside an unexpected shift's parentheses
are those between the nearest sequence positions on either side that
matched fragment ions of the spectrum. A single residue in parentheses,
`(A)[-0.0639]`, means the shift is pinned to that residue; a long range
means the spectrum has no fragment ions inside it.

## 3. Which modifications are shown by name

A modification is written by name when the search knew it in advance, and
as a mass when it did not:

| Kind | Given by | In the string | Columns |
|---|---|---|---|
| Fixed modification | TopPIC/TopMG `-f` (`C57`, `C58` or a modification file) | `(C)[Carbamidomethyl]` on every modified residue | `Fixed PTMs` |
| Protein N-terminal form | `-n` (`NME`, `NME_ACETYLATION`, `M_ACETYLATION`) | Acetylation as the prefix `[Acetyl]-`; methionine excision is not annotated, the proteoform simply starts at residue 2 | `Protein N-terminal form` |
| Variable modification | TopPIC `-b`, TopMG `-i` (a modification file) | `(S)[Phospho]` on the modified residue | `# Variable PTMs`, `Variable PTMs` |
| Truncation N-terminal modification | TopPIC `-b`, Position `N-term` (a variable PTM of the proteoform N-terminus, on a proteoform starting at residue 3 of the protein or later without a protein N-terminal form) | The prefix `[Gln->pyro-Glu]-QAEV...`, before any annotation of the first residue: `[Acetyl]-(K)[Acetyl]...` has an N-terminal and a side-chain acetylation on K | `# Variable PTMs`, `Variable PTMs` (`Acetyl:[N-term]`, first) |
| Unexpected mass shift | TopPIC `-s`, `-m`, `-M`; TopMG `-s` | `(...)[+79.9663]`, a signed mass over a residue range | `# Unexpected modifications`, `Unexpected modifications` |
| Characterized unexpected shift | TopPIC `-B` (a common modification file) with PTM characterization | The shift is replaced by the matching modification name, `(S)[Phospho]`, over the residues the characterization allows | `Variable PTMs`, `MIScore` |

The `Fixed PTMs`, `Unexpected modifications` and `Variable PTMs` columns
list the same annotations as `name:[position]` or `name:[first-last]`
(1-based positions in the proteoform), separated by `;`; for unexpected
shifts the "name" is the mass, for example `-0.0467:[7-33]`. A truncation
N-terminal modification is on the N-terminus rather than on a residue, so
its position is `N-term`: `Acetyl:[N-term];Acetyl:[1]` for
`[Acetyl]-(K)[Acetyl]...`.

**Where the names and masses come from.** The name in brackets is the
modification's abbreviation, and its mass is what TopPIC/TopMG add to the
residue mass:

- `C57` and `C58` are the built-in modifications `Carbamidomethyl`
  (57.021464 Da) and `Carboxymethyl` (58.005479 Da) on cysteine, and
  N-terminal acetylation is the built-in `Acetyl` (42.010565 Da). Their
  definitions are in the resource files `res/base_data/ptm_base.xml`,
  `mod_base.xml` and `prot_mod_base.xml`.
- A modification file (`-f <file>`, `-b`, `-B`, `-i`) supplies a name and a
  mass per line (`Name,Mass,Residues,Position,UnimodID`, see the TopPIC
  manual, section 5). The file is **not** the only source of the values,
  though: when the name is one of the abbreviations in
  `res/base_data/ptm_base.xml`, which are every Unimod (PSI-MS) name
  (`Phospho`, `Acetyl`, `Oxidation`, `Methyl`, `Carbamidomethyl`,
  `Cation:Na`, ...), the built-in entry is used and the mass in the file
  is ignored; when the two masses differ by more than 0.01 Da, TopPIC and
  TopMG print a WARNING saying so. Only a name that is not in that list is
  created from the file with the file's mass.
  So a line `Phospho,79.97,STY,any,21` in `common_mods.txt` still uses
  79.966331 Da, and a differently spelled name, say `Phosphorylation`,
  would use the mass from the file. To give a modification a different
  mass, use a name that is not in `ptm_base.xml`.
- `ptm_base.xml` holds every Unimod modification (1,560 in the Unimod
  release of 17 February 2026): the abbreviation is the Unimod name, the
  name the Unimod full name and the mass the sum of the AME2020 nuclide
  masses over the Unimod composition. Its only other entry is the empty
  `No PTM`. The names of earlier
  versions (`Carbamidomethylation`, `Oxoalanine`, `Allysine`,
  `Thiocarboxy`, `Cholinephosphoryl`, `TPQ`, and `carboxylation` and
  `Citrullination`, which are Unimod's `Carboxy` and `Deamidated`) are
  still accepted in modification files, FASTA sequences and earlier result
  files, and stand for the entry: the results show its present name.

## 4. The mass in brackets

An unexpected shift is written as a mass in Da with four decimals and an
explicit sign (`+` for positive shifts). It is the part of the observed
proteoform mass that the sequence with all its **named** modifications does
not explain: fixed, N-terminal and variable modifications are already
included in the residue masses, so the number never contains them. In the
example, the -116.0905 Da is in addition to the two carbamidomethylations
inside the parentheses.

Shifts close to zero (`-0.0467`, `+0.9477`) are usually not modifications
but the precursor mass error or a one-Dalton isotope error of the
deconvolution. Shifts close to a known modification mass (`+79.9663`,
`+42.0106`, `+15.9949`) can be turned into names by running TopPIC with a
common modification file (`-B`); the `MIScore` column then gives the
confidence of the assigned position.

The `Proteoform mass` column is the mass of the annotated proteoform,
including all named modifications and the mass shifts, and
`Adjusted precursor mass` is the precursor mass it was matched to.
Neither the string nor the two masses include an N-terminal label (the
`N-terminal label mass` column, 0 without a label): `Precursor mass` is
about `Adjusted precursor mass` plus the label.

## 5. Examples

The first two strings are from TopPIC result tables; the others are
constructed to show N-terminal acetylation and fixed and variable
modifications (`Previous amino acid` | `Proteoform` | `Next amino acid`):

```text
- | MKRTFQ(PSVLKRNRSHGFRARMATKNGRQVLAR)[-0.0467]RRAKGRARLTVSK | -
```
A complete protein (flanked by `-` on both sides) with a -0.0467 Da shift
located somewhere on residues 7 to 33 (`Unexpected modifications` =
`-0.0467:[7-33]`): a precursor mass error, not a modification.

```text
M | SLSTEATAKIVSEFGRDANDTGSTDVQVA(L)[+97.9219]LTAQINHLQGHF... | -
```
The initiator methionine is removed (`Protein N-terminal form` = `NME`,
`First residue` = 2, previous amino acid `M`), and a +97.9219 Da shift is
pinned to the leucine at position 30 of the proteoform.

```text
M | [Acetyl]-SLSTEATAKIVSEFGRDANDTGSTDVQVALLTAQINHLQGHF... | -
```
Methionine excision followed by N-terminal acetylation
(`Protein N-terminal form` = `NME_ACETYLATION`).

```text
K | (C)[Carbamidomethyl]SAEPKHKQ(S)[Phospho]RQG | -
```
A C-terminal fragment of a protein searched with `-f C57` and a variable
modification file containing `Phospho`: the cysteine is carbamidomethylated
(`Fixed PTMs` = `Carbamidomethyl:[1]`) and the serine phosphorylated
(`Variable PTMs` = `Phospho:[10]`).

```text
MKVRASVK(KL(C)[Carbamidomethyl]RN(C)[Carbamidomethyl]K)[-116.0905]IVKRDGVIRVI(C)[Carbamidomethyl]SAEPKHKQRQG
```
Three carbamidomethylated cysteines (fixed, `-f C57`) and one unexpected
shift of -116.0905 Da somewhere on `KLCRNCK`; the shift's parentheses
enclose two of the fixed modifications.

## 6. The ProForma column

The `ProForma` column of the same tables (and the `proforma` column of the
SQLite tables) writes the proteoform in ProForma 2.1, the notation of the
HUPO Proteomics Standards Initiative
(https://github.com/HUPO-PSI/ProForma), which other software can read.
The string is built from the proteoform itself
(`proteoform_util::getProFormaStr`), not from the `Proteoform` string:

- The residues are the proteoform's, from its first to its last residue,
  without the flanking residues.
- A PTM with a Unimod id (the `ptm_base.xml` entry, or the UnimodID column
  of a modification file) is written as its Unimod name and its Unimod
  accession, joined by `|` as synonyms: `[Oxidation|UNIMOD:35]`,
  `[Cation:Na|UNIMOD:30]`. The name is the abbreviation of the
  `ptm_base.xml` entry with that id, which is the Unimod name; a
  modification file's own name is not checked and not used for it (a line
  `Sodium,21.981943,DE,any,30` gives `[Cation:Na|UNIMOD:30]`). An id that
  `ptm_base.xml` does not have (one newer than its Unimod release) is
  written as the accession alone, `[UNIMOD:2200]`.
- A PTM without a Unimod id (a modification file's line with UnimodID -1,
  e.g. `MyPTM,12.345678,K,any,-1`) is written as its monoisotopic mass
  with its name as information, `[+12.3457|INFO:MyPTM]` (ProForma's
  `INFO` tag, which readers do not interpret). A name with `|`, `#` or unpaired square
  brackets is left out (`Cation:Mg[II]` is written).
- The protein N-terminal acetylations (`M_ACETYLATION`, `NME_ACETYLATION`)
  and a truncation N-terminal modification are N-terminal modifications,
  `[Acetyl|UNIMOD:1]-SEQ`.
- Fixed PTMs (`-f`), the PTMs written in the FASTA sequence and variable
  PTMs localized to one residue follow that residue:
  `C[Carbamidomethyl|UNIMOD:4]`, `M[Oxidation|UNIMOD:35]`.
- An unexpected shift is its mass with 4 decimals and a sign, on its
  residue (`L[+97.9219]`) or on its range of residues in parentheses
  (`(TLRL...NIM)[-236.9788]`); a variable PTM not localized to one residue
  is written on its range in the same way, with all the PTMs of the range,
  `(MAM)[Oxidation|UNIMOD:35][Oxidation|UNIMOD:35]`. Modifications of
  single residues can stand inside a range:
  `(MKAYC[Carbamidomethyl|UNIMOD:4]ERQGLSM)[Oxidation|UNIMOD:35]`.
- A PTM characterized by MIScore (`-B`) on several sites is written at its
  best site, with the score of every site with a positive score (ProForma's
  localisation scores); the label follows the first part of the tag:
  `K[Acetyl#s1(0.850)|UNIMOD:1]...K[#s1(0.150)]`.
- ProForma does not allow overlapping ranges: a range that would overlap
  one written before (from the left) is written as a set of positions, its
  first residue preferred, `P[+12.3456#g1]E[#g1]P[#g1]`, and a
  modification without a residue as one of unknown position,
  `[+12.3456]?SEQ`.

Everything except the joined name and accession, the localisation
scores, the sets of positions, the ranges and the unknown positions is
ProForma's base level (Level 1); those five belong to its Level 2 (a name
and an accession joined by `|` is its "joint representation", section 7.5;
a mass with an `INFO` tag is Level 1). A string with a PTM of Unimod is
therefore Level 2; `[Oxidation]` alone would be Level 1, but the
accession keeps the string exact if Unimod renames a modification.

A reader computes the mass of the proteoform from the string. For TopPIC
it is the `Proteoform mass` column, up to the rounding of the written
masses. TopMG's variable modifications keep the mass of the aligned gap
(an oxidation over a range may weigh 16.28 Da), which its `Proteoform
mass` column includes, while the ProForma string names the
modifications; its mass is then the theoretical one. Examples (`Proteoform` | `ProForma`):

```text
[Acetyl]-SGREGGKKKPLKQPKKQAKE(M)[Oxidation]DEEDKAFKQKQ...
[Acetyl|UNIMOD:1]-SGREGGKKKPLKQPKKQAKEM[Oxidation|UNIMOD:35]DEEDKAFKQKQ...

([Acetyl]-MPKGGR)[+365.1300]KGGH(KGR)[+167.0338]ARQYTSPEEIDAQ...
[Acetyl|UNIMOD:1]-(MPKGGR)[+365.1300]KGGH(KGR)[+167.0338]ARQYTSPEEIDAQ...

...KL(MKAY(C)[Carbamidomethyl]ERQGLSM)[Oxidation]RQIRFRFDGQ...
...KL(MKAYC[Carbamidomethyl|UNIMOD:4]ERQGLSM)[Oxidation|UNIMOD:35]RQIRFRFDGQ...
```

