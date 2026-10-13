# Project conventions and migration notes

Guidance for working on this codebase: the code-wide conventions, the
dependencies, the build and install rules, the manuals and the source
layout. The notes on each tool are in `.claude/tools/` (next section).

## Tool notes (`.claude/tools/`)

The notes on each tool's algorithms, options and results are in one file
per tool under `.claude/tools/`. They are not loaded automatically: read
the file of the tool you work on before changing its code, and update it
in the same commit as the code (as this file was before 2026-10-12).

- `.claude/tools/topfd.md` — TopFD and TopDIA:
  - Quantized ONNX models
  - Envelope refinement and the core-peak intensity fit
  - Deconvolution methods
  - FLASHDeconv MS/MS deconvolution
  - Greedy MS1 deconvolution
  - Greedy MS/MS deconvolution
  - Precursor isotope correction from fragment pairs
  - PTCR spectra
  - MS3 spectral deconvolution
  - Activation of MS/MS spectra
  - Isotope distributions
  - The MS1 peaks of the database and the 3D levels
- `.claude/tools/toppic_search.md` — TopPIC tolerances, filters and searches (also the index files topindex writes):
  - Precursor and fragment mass tolerances
  - Precursor error tolerance in Dalton
  - Precursor +-1 Da error per identification stage
  - Precursor refinement of the variable PTM search
  - Precursor refinement of the shift searches
  - The filters' mass list
  - Tag filter
  - Second alignment of the shift searches
  - Large-shift penalty of the shift alignments
  - Range of a variable PTM shift
  - Truncation N-terminal modifications
  - Index files version 2, the combined filter pass and the E-value stage (toppic_test experiment 150)
- `.claude/tools/toppic_evalue.md` — TopPIC E-values and FDRs:
  - E-values of PrSMs found at a +-1 isotope offset
  - E-values of the variable PTM PrSMs with modified residues
  - Tolerance of the complementary PRM peaks in the E-values
  - Deamidation or the -1 isotope offset
  - Rescoring of the PrSMs
- `.claude/tools/toppic_output.md` — TopPIC after the search (also topcalib, and the PrSM tables topmg writes too):
  - Post mass matching
  - ProForma strings
  - Average masses in msalign files
  - Proteoform classification
  - PTM localization by a variable PTM search
  - The second precursor search
  - Two-pass m/z recalibration
- `.claude/tools/topmg.md` — TopMG:
  - Graph alignment tolerance
  - Protein N-terminal forms in TopMG
  - PTMs written in FASTA sequences in TopMG
  - Precursor refinement of TopMG's graph alignment
- `.claude/tools/topdiff.md` — TopDiff:
  - TopDiff
- `.claude/tools/toplib.md` — toplib_build and toplib_search:
  - toplib_build
  - toplib_search

## Removed functions in `str_util.hpp` (use the standard library)

These thin wrappers were removed because they merely duplicated the C++ standard
library. Do **not** re-add them; call the standard function directly.

| Removed | Use instead |
|---|---|
| `str_util::toString(int)` | `std::to_string(i)` |
| `str_util::toString(size_t)` | `std::to_string(n)` |
| `str_util::scientificToDouble(s)` | `std::stod(s)` (also parses scientific notation, e.g. `"1.5e-3"`) |

Still present (kept on purpose): `toString(bool)` ("true"/"false"),
`toString(double)` and the `*ToString` formatters (custom precision logic),
`trim`, `split`, `rmComment`, and `endsWith`. `endsWith` is kept rather than
using C++20 `std::string::ends_with`, because the project targets **C++17**.

Note: `str_util` no longer depends on Boost — `trim`/`split`/`rmComment` are
implemented with the standard library; `split` reproduces the old
`boost::split + is_any_of` semantics (empty fields preserved, `""` -> `{""}`).

## Removed functions in `file_util.hpp` (use `std::filesystem`)

Only the two exact one-call `std::filesystem` duplicates were removed:

| Removed | Use instead |
|---|---|
| `file_util::exists(p)` | `std::filesystem::exists(p)` |
| `file_util::rename(a, b)` | `std::filesystem::rename(a, b)` |

The other `std::filesystem`-backed helpers are **kept** as a `std::string`
facade (they centralize error handling/logging and keep callers off `fs::path`):
`directory`, `absoluteName`, `absoluteDir`, `filenameFromEntirePath`,
`basenameFromEntirePath`, `basename`, `copyFile`, `delDir`, etc.

Watch out: `file_util::basename(s)` is **not** `fs::path::stem()` — it keeps the
directory and strips only the last extension (`/a/b/c.txt` -> `/a/b/c`), i.e. it
is `fs::path(s).replace_extension("")`. Use it, don't replace it with `stem()`.

## Migrating XML code from Xerces-C++ to pugixml — possible bites

`common/xml` is fully on **pugixml** (Xerces is gone). When migrating code that
still uses Xerces (`appendXml`/`toXmlElement`/`parseXml`, writers, base
`initBase`), watch for these:

- **Element handles, not pointers.** `XmlDOMElement` is now `pugi::xml_node` (a
  value handle), not `xercesc::DOMElement*`. Pass by value; test absence with
  `if (!node)`, not `== nullptr`.
- **No detached nodes.** `XmlDOMDocument::createElement(tag)` (detached, attach
  later) is gone. Use `addElement(parent, tag)` which creates the child **already
  attached** to `parent` and returns it. Build the subtree under the returned
  node. `createTextNode` folds into `addElement(parent, tag, value)`.
- **No serializer object.** Use `xml_dom_util::writeToString(node)` (a node
  method underneath), not `writeToString(serializer, node)`.
- **No `DOMImplementation` singleton.** `XmlDOMImplFactory` / `createSerializer`
  are gone; `xml_dom_impl::createDoc(root)` returns a `unique_ptr<xml_document>`.
- **Parsing.** `XmlDOMParser::parse(file)` / `parseStr(buffer)` return
  `std::unique_ptr<pugi::xml_document>`; there is no `MemBufInputSource` and no
  `XMLPlatformUtils` init. Construct `XmlDOMDocument` from the returned doc.
- **RAII.** `XmlDOMDocument` owns its document; there is no `release()`.
- **Errors.** pugixml has no SAX error handler. Parse errors come back as an
  `xml_parse_result`; `xml_dom_err_handler::checkParseResult(result, name, src)`
  logs and throws (reporting line:column when `src` is provided).

### The semantic gotcha (most important)

Xerces `getElementsByTagName(tag)` searches **all descendants** (any depth, in
document order). pugixml `node.children(tag)` returns only **direct children**.
The migrated `xml_dom_util::getChildElement/getChildValue/getChildCount` use
**direct children**. An audit confirmed this is safe for TopPIC's existing XML
(it is read level-by-level and every collection is wrapped in a dedicated
`*_list` element), so direct-children and descendant search coincide there.

But if you write/query XML where the **same tag nests at multiple depths under
the queried element**, the two differ (count and indexing), e.g. counting
`<name>` on a `<mod>` whose subtree also has `<name>` deeper. If you genuinely
need the descendant search, use pugixml XPath: `parent.select_nodes(".//tag")`
(sort for document order). Also note `getTextContent` (concatenates all
descendant text) became `text().as_string()` (first PCDATA) — identical for leaf
value elements, different for mixed content. pugixml assumes UTF-8 and has
minimal XML-namespace support.

## Header include order (required)

All `.cpp` and `.hpp` files use the standard (Google-style) order. Each group is
**alphabetical** and separated from the next by **one blank line**:

1. **The related header** — for a `.cpp`, its own `.hpp`, on the first line.
2. **C++ standard library** headers (`<string>`, `<cmath>`, `<cstddef>`, ...).
3. **Platform / system headers**, including `#if`-guarded blocks (`<windows.h>`,
   `<unistd.h>`, `<mach/mach.h>`, ...).
4. **Third-party** headers (`<pugixml.hpp>`, `<boost/...>`).
5. **Project** headers (`"common/..."`).

A `.hpp` follows the same order without step 1. Putting the related header first
keeps each header self-contained. Include hygiene is checked with
include-what-you-use (configure with `-DTOPPIC_ENABLE_IWYU=ON`, or run
`iwyu_tool -p <build>`); keep new files IWYU-clean.

## Parameter passing (smart pointers and handles)

The `...Ptr` typedefs are `std::shared_ptr` aliases. A `shared_ptr` is 16 bytes
and a by-value pass does an **atomic** refcount inc/dec on every call, so:

- **Read-only `shared_ptr` params → `const&`.** `void f(const ModPtr &p)`, not
  `void f(ModPtr p)`. Same for `std::shared_ptr<T>` written out. This is applied
  throughout `seq` and `common/base`; keep new code consistent.
- **Sink params that store the pointer → by value + `std::move`.** A setter or
  constructor that *keeps* its argument takes it by value and moves it into the
  member (`void setX(ModPtr p) { x_ = std::move(p); }`). Do **not** "fix" these
  to `const&` — that would force a copy on assignment.
- **`...PtrVec` (and other container) params → `const&`** when read-only.
- **`XmlDOMElement` stays by value** — see the pugixml section above. It is a
  trivially-copyable 8-byte handle (`pugi::xml_node` wraps a single pointer), so
  `const&` would only add a layer of indirection. The same reasoning applies to
  any thin handle/iterator-like type.

Scalars (`int`, `double`, `bool`, ...) stay by value. The one risk `const&`
introduces and the compiler will **not** catch: if a function mutates the
container its `shared_ptr` argument was passed from while still using the
argument, the reference can dangle — pass such an argument by value.

- **Getters of `shared_ptr` members on hot paths → `const&` too**, and never
  copy a `shared_ptr` to a *shared* object (base-table residues, amino
  acids, PTMs, `BasePeakType::ORIGINAL`, `ProteoformType::*`) inside a
  per-residue or per-peak loop. The atomic refcount update lands on a cache
  line every search thread hammers, and it collapsed the threaded zero/one
  shift searches to a 2x speedup on 16 threads: per-residue copies in
  `residue_util::convertStrToResiduePtrVec` (a temporary `Residue` per
  residue in `ResidueBase::getBaseResiduePtr`), `ResidueSeq::getResiduePtr`
  and `Residue::getAminoAcidPtr/getPtmPtr` (by value, hit by `BpSpec`), and
  `PrmPeak::getBaseTypePtr` (copied per iteration of
  `diag_pair_util::compDiagPair`). Fixing those gave 6.8x and 11x. Compare
  such pointers with `.get()` or `==` on references; `AminoAcidBase::
  getAminoAcidPtrByOneLetter`, `PtmBase::getPtmPtrByAbbrName` and
  `getEmptyPtmPtr` return `const&` for this reason.

## Vendored htslib (`seq` indexed-FASTA access)

`seq/fasta_index_reader` uses htslib's `faidx` API for random access into
`.fai`-indexed FASTA files. A **trimmed** copy of htslib (only `faidx`, `bgzf`
and `hfile` — the C sources plus their headers) is vendored under `ext/htslib`,
mirroring how the upstream TopPIC tree carries it. It is built by `CMakeLists.txt`
as a small static library `htslib` (links `ZLIB::ZLIB` and `Threads::Threads`,
compiled with `-w` and PIC) and linked **PRIVATE** into `toppic_common`; the
`ext/` directory is a **PUBLIC** include root because the vendored sources
include themselves as `"htslib/<name>.h"` and `fasta_index_reader.hpp` exposes
`<htslib/faidx.h>`.

Notes:
- It is third-party C code: do not reformat it or hold it to this project's
  include-order / IWYU rules. The `-w` flag deliberately silences its warnings.
- do not add htslib calls outside `seq` without reason.

## Other external dependencies

Besides pugixml (XML) and the vendored htslib (indexed FASTA), the layers below
pull in a few system libraries, all found via `find_package`/the default include
path and linked into `toppic_common`:

- **SQLite3** (`find_package(SQLite3)` → `SQLite::SQLite3`) — used by
  `sql/sql_util` and `ms/mzml/mzml_ms_sql_writer`. The latter is the single-file
  store for deconvoluted MS1/MS2 spectra (it replaced an earlier one-JSON-file-
  per-scan writer): it inserts via reused prepared statements batched into
  chunked transactions under bulk-load PRAGMAs, so it scales to many scans.
  **Interface note:** `mzml_ms_sql_writer` is the stateful class
  `MzmlMsSqlWriter` (construct once, call `writeMs1`/`writeMs2` per spectrum, let
  it destruct to flush the final transaction), **not** the old
  `mzml_ms_sql_writer::writeMs1/writeMs2(sqlite3*, ...)` free functions — the
  batching needs state a per-call free function can't hold. The caller still owns
  the connection (opens it, creates the schema, closes it). `ms1_peak` is a
  `WITHOUT ROWID` table keyed on `(LEVEL, spec_id, peak_id)` (toppic_test
  experiment 165): the writer gives a peak `LEVEL` 0 (-1 for intensity <= 0),
  and with topfd `--sql-3d` it writes no `ms1_peak` rows (its constructor's
  `write_ms1_peaks`), because `MsSqlWriter` writes them with the 3D levels. When migrating the
  upstream topfd driver that wrote spectra, adopt the object; do not reintroduce
  per-call free-function wrappers (they would restore the per-scan
  transaction/prepare cost this rewrite removed).
- **Boost** — `find_package(Boost ... COMPONENTS filesystem iostreams thread
  chrono serialization program_options)`. `ms/mzml`'s pwiz reader needs
  `filesystem`/`iostreams`/`thread`/`chrono`; `filter/massmatch/
  mass_match` needs `serialization` (its binary index archive — linked PUBLIC
  into `toppic_common` because the header exposes `boost::serialization::access`);
  the GUI argument parsers need `program_options`. `search/graph` additionally
  uses the **header-only** Boost Graph Library (no link). Prefer `std::mutex` etc.
  over `boost::*` in our own code (the filter/index mutexes were converted from
  `boost::mutex`; only `mass_match`'s serialization genuinely needs Boost).
  **Do not add `system` back** to the component list or `Boost::system` to a
  link line: Boost.System is header-only since 1.69, and Boost >= 1.90
  (Ubuntu 26.04) ships neither the stub library nor a `boost_system` CMake
  config, so `find_package` fails with "Could not find ... boost_system";
  Ubuntu 26.04 has no `libboost-system-dev` package either. (Boost
  headers pulled in through the third-party `ext/` tree are silenced by the
  system-include treatment below; if you include a Boost header directly in our
  own code and it warns — e.g. uBLAS still deriving from the C++17-deprecated
  `std::iterator` — wrap it in `#pragma GCC diagnostic ignored
  "-Wdeprecated-declarations"`.)
- **Compiler default** — `CMakeLists.txt` defaults to clang/clang++ when found,
  except on Windows hosts (MSYS2), where it leaves CMake's `cc`/`c++` default
  so the compiler matches the one the MSYS2 packages were built with (GCC in
  UCRT64/MINGW64, clang in CLANG64). Mixing clang objects with GCC-built
  static Boost fails at link time with "duplicate section ... has different
  size" warnings and "multiple definition" errors for COMDAT variables that
  clang places in `.bss` and GCC in `.data`. Keep that `CMAKE_HOST_WIN32`
  guard; `docs/compile/windows_build.md` tells users to install `...-gcc`, not `-clang`.
- **Qt6** — `find_package(Qt6 COMPONENTS Widgets Core Gui)` backs the `src/gui`
  desktop executables (see the source-layout note). Only the GUI targets use it,
  via per-target `AUTOMOC`/`AUTOUIC`/`AUTORCC`; the `toppic_common` library has
  no Qt dependency.
- **ProteoWizard (pwiz)** — version 3.0.26278, vendored under `ext/pwiz` with
  only the files `ms/mzml/pw_ms_reader` needs (the sources the linker takes and
  the headers they reach; `ext/README.md` lists the rule, the local changes and
  how to update). Every `.cpp` under `ext/pwiz` is compiled into a static
  `pwiz` library against Boost, zlib and zstd (pwiz's binary data encoder
  decodes zstd-compressed arrays) with `WITHOUT_MZ5` and `WITHOUT_MZMLB` (so no
  HDF5; mz5 and mzMLb are not supported). `ms/mzml/pw_ms_reader` uses it to
  read mzML/mzXML. pwiz uses its own `minimxml` parser, **not** Xerces. It needs
  a few headers of pwiz's `boost_aux` that Boost itself does not have (the
  `boost::enums` library, the singleton utility, `boost/nowide/system.hpp`);
  those are vendored under `ext/boost`, and everything else, Boost.Nowide
  included, comes from the system Boost. The nowide parts pwiz uses are
  header-only off Windows; on Windows `pwiz` links `Boost::nowide` (the UTF-8
  console streams `boost::nowide::cerr` etc. of pwiz's `Stream.hpp`, and the
  file buffer), which `find_package(Boost)` requests only there. Treat `ext/pwiz` and
  `ext/boost` as third-party: do not reformat them or hold them to the
  IWYU/include-order rules.
  Their `.cpp` files are built with `-w`; the `ext/` directory is added to every
  target as a **`SYSTEM`** include (`-isystem`), so warnings from vendored
  headers are silenced even when our own `-Wall` sources include them.
- **ONNX Runtime** — `topfd/envcnn` and `topfd/ecscore/score` run the EnvCNN /
  ECScore neural-network models through the ONNX Runtime C++ API. A prebuilt
  `libonnxruntime.so.1.23.2` (SONAME symlink `libonnxruntime.so.1`; glibc >= 2.27; **stored with Git LFS**, rule
  `ext/onnx/libonnxruntime.so.*` in `.gitattributes`, with a `.gitignore`
  exception to the `*.so.*` rule) and the C/C++ API headers are vendored under
  `ext/onnx` (sources include it as `"onnx/onnxruntime_cxx_api.h"`, covered by
  the `SYSTEM` `ext/` include); it is linked as an `IMPORTED` shared object,
  `PRIVATE`, into `toppic_common`. The vendored `.so` is Linux x86-64 only;
  on macOS and Windows `CMakeLists.txt` instead locates the system package
  (Homebrew's `libonnxruntime.dylib`; MSYS2 UCRT64's `libonnxruntime.dll.a`)
  with `find_library` (overridable with `-DONNXRUNTIME_LIBRARY=...`) as an
  `UNKNOWN IMPORTED` target, and still compiles against the vendored 1.23
  headers. Windows also links `ws2_32` into the vendored htslib and pwiz
  (Winsock), builds the GUI targets as `WIN32_EXECUTABLE`, and installs `res/`
  next to the executables because `getResourceDir()` does not consult the
  compiled-in shared dir there. Treat it as third-party, with one local
  patch: `onnxruntime_c_api.h` defines `ORT_API_CALL` as `__stdcall` on
  Windows (upstream writes `_stdcall`, which GCC rejects under the strict
  `-std=c++17` that `CMAKE_CXX_EXTENSIONS OFF` selects; MSYS2 patches its
  packaged header the same way). Re-apply it when re-vendoring the headers.

`file_util::getResourceDir()` looks for the resources in `<exec_dir>/res`
(installed layout on macOS/Windows), then `<exec_dir>/../res` (the build tree:
executables in `<repo>/bin`, resources in `<repo>/res`, so nothing needs to be
copied or symlinked to run from `bin/`), then the compiled-in shared dir
(installed layout on Linux). Keep that order if you touch it.

## Install and uninstall (`make install` / `make uninstall`)

`make install` is defined by three rules in `CMakeLists.txt`, all relative to
`CMAKE_INSTALL_PREFIX` (default `/usr/local`), so `DESTDIR` staging and a
custom prefix both work:

- `install(TARGETS ${TOPPIC_EXECUTABLES} ...)` puts the nine CLI tools and
  the six `*_gui` counterparts (toplib_build, toplib_search and topcalib have
  none) in `<prefix>/bin`. **When you add a new executable
  target, append it to the `TOPPIC_EXECUTABLES` list** (end of the file) so it
  is installed and gets the rpath below — nothing else is needed. The one
  exception is a developer tool outside the default target, such as
  `theo_patt_gen` (`EXCLUDE_FROM_ALL`), which stays out of the list.
- On Linux the vendored `ext/onnx/libonnxruntime.so.1.23.2` and its SONAME symlink (which every
  executable loads at run time, via `toppic_common`) is installed to
  `<prefix>/${CMAKE_INSTALL_LIBDIR}` (= `lib/toppic`), and the executables get
  `INSTALL_RPATH "$ORIGIN/../lib/toppic"` so the installed tree is relocatable.
  In the build tree CMake already points the rpath at `ext/onnx`, so `bin/`
  works without installing. Keep the `CMAKE_INSTALL_LIBDIR` value and this
  rpath in sync.
- `install(DIRECTORY res/ ...)` copies the runtime resources into
  `<prefix>/${CMAKE_INSTALL_DATADIR}` (= `share/toppic`), which is also the
  compiled-in `TOPPIC_SHARED_DIR` the binaries look in.

`make uninstall` is a custom target (CMake has no built-in one). It runs
`cmake -P <build>/cmake_uninstall.cmake`, configured from
`cmake/cmake_uninstall.cmake.in`, which reads the `install_manifest.txt` that
`make install` writes into the **build directory**, deletes every file listed
there (honouring `DESTDIR`) and then removes any directory under the prefix
that was left empty (the manifest lists files only). Consequences:

- It must be run from the same build directory as the install; without a
  manifest it fails with a clear error rather than guessing.
- Anything added through a normal `install(...)` rule is uninstalled
  automatically — do not add per-file removal code to the script.
- Running it twice is harmless (missing files are reported, not errors).

## Platform build documentation

The build instructions are split by platform. `README.md` has one "Building
TopPIC Suite" section with a link per platform; the details live in
`docs/compile/`, each
starting with a "System requirements" list for that platform (compiler,
CMake, Boost components, libraries, where ONNX Runtime comes from, Git
LFS) followed by the same six numbered steps:

- `docs/compile/ubuntu_build.md` — the primary platform: `apt-get` package
  list (also `scripts/install_ubuntu_packages.sh`), clone with Git LFS,
  configure/build, run from `bin/`, install, uninstall.
- `docs/compile/windows_build.md` — MSYS2 **UCRT64** shell, Ninja generator, GCC (not
  clang; see the compiler-default note above), `pacman` package list, system
  ONNX Runtime, `res/` installed next to the executables.
- `docs/compile/redhat_build.md` — RHEL 10 (and Rocky/AlmaLinux): CodeReady Builder +
  EPEL for `pugixml-devel`, `dnf` package list, otherwise identical to Ubuntu.
- `docs/compile/macos_build.md` — Xcode Command Line Tools + Homebrew, Homebrew
  `onnxruntime` (nothing is vendored for macOS; `-DONNXRUNTIME_LIBRARY`
  override), Qt6 hint, resources installed to `<prefix>/bin/res`.

When you change anything platform-specific in `CMakeLists.txt` — a required
package, the compiler choice, how a library is located, the install layout —
update the matching doc in the same commit; the
docs are hand-written and nothing checks them against the build.

## User manuals (`docs/manual/`)

One manual per command-line tool, linked from README.md's "Manuals" section:
`topfd_manual.md`, `topdia_manual.md`, `topindex_manual.md`,
`toppic_manual.md`, `topmg_manual.md`, `topdiff_manual.md`,
`toplib_build_manual.md`, `toplib_search_manual.md`,
`topcalib_manual.md`. Each follows the
same layout: overview and usage line, input (including the companion files a
tool needs, e.g. TopPIC's `_ms2.feature` and `.sqlite`), the pipeline steps
as the console prints them, the output files with their columns, option
tables of the **visible** options with defaults (hidden options — the ones
the parser accepts but `-h` does not print — are deliberately not
documented in the manuals), examples. They were written from the
`console/<tool>_argument.cpp` parsers, the `*_process.cpp` drivers and the
writers, and nothing checks them against the code, so:

- When you add, rename, hide/unhide or change the default of an option in a
  `<tool>_argument.cpp`, or change an output file name, format or the
  columns of a writer, update that tool's manual in the same commit. The
  option tables mirror the `-h` text; keep the two saying the same thing.
- Do not document an option the parser does not accept. Cross-tool facts
  that several manuals repeat and must stay consistent: result-file naming
  (`_post_ms2` after TopPIC's default post mass matching, plain `_ms2` for
  TopMG and with `-D`), which tool needs which TopFD file, the modification
  file format (`Name,Mass,Residues,Position,UnimodID`), and the shared PrSM
  table columns (documented once in the TopPIC manual, referenced by TopMG).
- `topfd_manual.md` specifics: its examples use no `-a` (activation
  defaults to `FILE`), and `-T`-mode claims must be checked against
  `topfd_single_process.cpp` (that mode ignores the MS1/feature options).
- TopDIA's `-y`/`-Y` options (`-p`/`-P` until 2026-10-08) set
  `TopfdPara::seed_env_inte_corr_tole_cutoff_`
  (copied into `EcscorePara` by `env_coll_detect::processMs1/processMs2`)
  before each feature-detection stage; both default to 0.5, the value the
  detector had hard-coded, so default runs are unchanged. TopPIC's
  `-N`/`--n-terminal-label` works (its handler used to test the wrong key); it
  is a visible option and a toppic_gui field ("N-terminal label mass (Da)",
  Advanced Parameters tab, `nTermLabelEdit`; `gui/util/command.cpp` emits
  `-N` only for a non-zero value). (TopDiff resolves result
  files itself: `<base>_post_ms2_<tool>_proteoform.xml` first, then
  `<base>_ms2_...`, for either `_ms2.msalign` or `_post_ms2.msalign` input —
  `merge/feature_sample_merge.cpp`; keep the manual's section 1 in step.)

`docs/manual/proteoform.md` (linked from the TopPIC and TopMG manuals and
the README) explains the annotated proteoform string of the `Proteoform`
column: it is built by `Proteoform::getProteoformMatchSeq` (`seq/
proteoform.cpp`, labels from `MassShift::getAnnoStr`), and modification
names/masses resolve through `PtmBase::getPtmPtr`, which prefers the
`res/base_data/ptm_base.xml` entry over a modification file's mass when the
abbreviation matches. Update it when that code or the table columns change.

`res/README.md` describes the runtime resource files (which class loads
each `base_data` file, the two ONNX models, the Git LFS-tracked files);
update it when a resource file is added, removed or repurposed.

## Source layout

`src` is the include root, so headers are included by their path from `src`
(`common/...`, `seq/...`, `ms/spec/...`, `para/...`, `sql/...`). All of the
following compile into the single `toppic_common` **OBJECT** library (the
`COMMON_SRCS` glob in `CMakeLists.txt` lists each directory) — its object files
are linked directly into each executable, so **no `libtoppic_common.so`/`.a` is
built**; the executables are self-contained. Each layer depends only on the ones
above it, never the reverse:

- `src/common` — foundation: `base`, `util`, `xml`, `thread`, `emass`
  (isotope distributions; see `.claude/tools/topfd.md`).
- `src/sql` — thin SQLite helper (`sql_util`: `execSql`, `prepareSql`,
  `stepAndReset`, `hasTable`, all abort on error — use them rather than
  re-implementing the prepare/step boilerplate), `sql_db_copy` (toppic's
  copy of TopFD's database, `<base>_toppic.sqlite`, plain or corrected by
  `--mass-calib`: backup API, toppic's tables dropped, written with no
  journal under a temporary name, closed in WAL mode (TopFD's journal mode:
  in rollback mode a reader of the copy, such as a viewer, made toppic's
  next write stop with "database is locked", toppic_test experiment 156)
  and renamed when complete; `removeDb` deletes a database with its
  -journal, -wal and -shm files) plus
  `ms_sql_writer`, the
  MS1-peak levels for 3D visualization: `MsSqlWriter` builds level 0 (every
  MS1 peak of intensity > 0) and the down-sampled levels 1..n, and writes
  them as the `LEVEL` column of `ms1_peak` (the highest level that keeps the
  peak; the levels are nested) plus one `CONFIG` row per level, on a
  caller-owned connection like the topfd writer below; readers compute the
  colour bucket (before toppic_test experiment 165 the levels were the
  tables `PEAKS0..n` with a `COLOR` column; see `.claude/tools/topfd.md`).
  topfd `--sql-3d` (in `deconv_ms1_process`) writes them into the topfd
  `.sqlite` after MS1 deconvolution.
- `src/para` — analysis parameters (`sp_para`, `peak_tolerance`, `prsm_para`).
- `src/seq` — sequence/proteoform layer.
- `src/ms` — mass-spectrum layer: `spec` (peaks/spectra/msalign), `msmap`,
  `factory`, `env` (envelope detection), `feature`, `mzml`, `calib` (the
  two-pass m/z recalibration of topcalib and toppic `--mass-calib`, with
  the corrected database of `sqlite_mass_calib`, see
  `.claude/tools/toppic_output.md`).
  `MsAlignReader` has a four-argument constructor taking an EnvCNN score
  cutoff (`SpPara::getEnvCnnCutoff()`, whose own default is 0.2, but both
  toppic and topmg set the key from their `-F` (`--filter-by-score`; the
  fourth msalign column is TopFD's EnvCNN score or TopDIA's pseudo score;
  `--filter-by-envcnn` until 2026-10-12; the parameter block, topmg_gui and
  toppic_gui (Advanced Parameters tab, `envCnnCutoffEdit`, since 2026-10-12)
  say "Confidence score cutoff", "EnvCNN score cutoff" before),
  which defaults to 0.4 in
  each): masses whose fourth
  msalign column is below it are skipped, keeping the file's peak ids. The
  filter, search and E-value/p-value stages (`filter/*`, `search/*`,
  `stat/tdgf`, `stat/mcmc`, `search/graph/spec_graph_reader`) construct their
  readers with it, so they work on the filtered masses; the PrSM output
  stages (`prsm_reader_util`, coverage, localization, tables) deliberately
  read the unfiltered spectra. Use the cutoff constructor in a new
  identification stage, and the three-argument one in a new output stage.
  Because the search stores matched mass/fragment numbers counted on the
  filtered masses, toppic and topmg run `prsm/prsm_match_num_recount`
  (writing `<tool>_recount`) to recount them against the full spectra --
  topmg on `topmg_raw_prsm` before its clustering, toppic on
  `toppic_cluster_prot` after its clustering and cluster protein correction (`prsm_match_num_recount::process` and
  `prsm_reader_util::readPrsmsWithSpectra` take an EnvCNN cutoff whose
  default, 0, keeps every mass; both drivers use the default); the number
  of masses in the tables is not stored in the PrSM but computed from the
  attached full spectra at output.
  Spectrum ids (`SPECTRUM_ID`) need not start at 0 or be consecutive: a
  subset of a file keeps the ids of its spectra, and `msalign_frac_merge`
  adds `SpPara::getMaxSpecNumPerFile()` (10,000,000) per fraction. Code that
  keeps one entry per spectrum indexes it by the spectrum's position in the
  file (`msalign_util::readSpIds`, the ids in file order), never by its id
  up to `getSpNum`: the filters' candidate lists (`SimplePrsmMemMerge`,
  `DiagFilterProcessor`) once allocated one list per id below the spectrum
  number and dropped the candidates of every larger id without a word, so
  TopPIC found nothing in a file whose ids start above 0 (issue #598,
  toppic_test experiment 100).
- `src/topfd` — the TopFD deconvolution/feature-detection layer, built on `ms`:
  `common` (`topfd_para` config + the `topfd_process`/`topfd_single_process`
  orchestrators), `dp` (dynamic-programming envelope assignment), `envcnn` and
  `ecscore/score` (the EnvCNN / ECScore neural-network scorers, via ONNX
  Runtime), `deconv`, and `ecscore` (`env`/`env_set`/`env_coll`/`para`/`score`).
  Note two folder-level cycles handled by the single-library glob: `dp <-> deconv`
  and `ecscore/env_coll <-> ecscore/score`. `deconv` constructs the
  `MzmlMsSqlWriter` (see the SQLite note) once per run, shared across its worker
  threads; it no longer writes per-scan JSON.
  TopFD has two entry points, chosen in `console/topfd.cpp` by
  `TopfdPara::isTextPeakList()` (`-T`): `topfd_process` deconvolutes
  mzML/mzXML files (MS1 deconvolution → feature detection → MS/MS
  deconvolution, one fraction per FAIMS voltage), and `topfd_single_process`
  deconvolutes one MS/MS spectrum from a text peak list (`m/z intensity` per
  line, separated by spaces or a tab; `readPeakFile` skips the lines that
  do not start with a digit, so an MGF file of one spectrum is read as
  one, and the peaks of all spectra of a longer MGF file are merged, with
  a warning when it has more than one `BEGIN IONS` line)
  with `DeconvSingleSp`, writing `<input>_ms2.msalign`,
  `<input>_ms2.env` and `<input>.sqlite`. **Both** name their outputs through
  `TopfdPara::setMzmlFileNameAndFaims` (base name = input path minus
  extension, `_<voltage>` appended for FAIMS; it also creates the SQLite
  database when enabled) — do not hand-roll output names or a second
  `createSqlDb` call in a driver. In topfd_gui the checkbox "Text peak list
  input" (`textPeakListCheckBox`, Additional settings, third column, row 3)
  sets `setTextPeakList` and `gui/util/command.cpp` emits `-T`; the Add
  button then offers `*.txt` and `*.mgf` or all files, and Start refuses a
  list holding an mzML/mzXML file (none of its lines starts with a digit,
  so `readPeakFile` would return an empty spectrum).
- `src/topdia` — the TopDIA pseudo-spectrum layer built on `topfd`/`ms`:
  `common` (`topdia_para` + `topdia_process`) and `pseudo_spec`
  (`mzrt_feature`, `pseudo_peak`, `pseudo_spectrum`, `generate_pseudo_spectrum`).
- `src/prsm` — proteoform-spectrum-match layer (`Prsm`/`SimplePrsm`, readers,
  writers, FDR, clustering, coverage). `prsm_reader_util::readPrsmsWithSpectra`
  reads a PrSM file and attaches the deconvoluted/refined spectra from the
  msalign file; both `prsm_match_table_writer` (the TSV tables) and
  `prsm_sql_writer` (the same columns as `prsm`/`proteoform`/`protein` tables, plus
  `prsm_mass_shift` and `prsm_protein_match`) are built on it.
  `prsm_second_prec` is the second precursor search and `prsm_mass_calib`
  reads the calibrants of `--mass-calib` (see
  `.claude/tools/toppic_output.md`).
  `prsm_sql_output::write` is the one entry point that writes those tables
  plus `seq/fasta_sql_writer`'s `fasta_seq` table into the database it is
  given when it exists (toppic: its `<base>_toppic.sqlite`; topmg: the
  topfd `.sqlite` of the spectrum file); toppic and topmg both call it after
  their protein table, passing their own cutoff file extensions — do not duplicate
  the open/write/close sequence in a driver. Keep the TSV and SQL columns in
  step when changing either. `prsm`/`simple_prsm`/`expected_value`/
  `peak_ion_pair` carry the pugixml XML (de)serialization (`toXmlElement`/`toXml`/
  `appendXml` take a parent `XmlDOMElement` and return the attached node).
  `search_fasta_match` (the "# Protein hits" / `prsm_protein_match` lookup
  of every protein containing a proteoform's sequence, run several times
  per PrSM by the table and SQL writers) is a 4-mer position index over
  the concatenated database: a query verifies only the positions of its
  rarest 4-mer, which gave the same results as the previous
  `std::string::find` scan of every protein and ran 59x faster on 1216
  real proteoforms against the 20k-protein human database (6.9 ms -> 0.12
  ms per query; queries shorter than 4 residues fall back to the scan).
  The index costs about 4 bytes per residue on top of the sequences.
- `src/filter` — proteoform filtering (`diag`/`index`/`massmatch`/`mng`/
  `oneptm`/`varptm`/`zeroptm`). `massmatch/mass_match` (de)serializes its binary
  index with Boost.Serialization. The file-global `serial_mutex` (and the other
  filter/index mutexes) live in **anonymous namespaces** — they previously had
  external linkage and the same name in several TUs, which only becomes a
  multiple-definition error once everything compiles into one library.
  `filter/index` is topindex's side: `zero_ptm_index` and `one_ptm_index`
  build, per database block, `MassMatch` objects from the same factory
  calls with only the manager's `max_proteoform_mass_` / `filter_scale_`
  differing (zero-shift 100000/100, one-shift 50000/100). Because the
  construction loops stop at the first mass beyond the column limit and the
  prefix/suffix masses are nondecreasing, a one-shift index is the
  zero-shift index of the same kind cut at the one-shift limit
  (`MassMatch::truncateColumns`, guarded by `canTruncateColumns`, which
  checks that the rows were nondecreasing). topindex therefore runs
  `joint_index::process`: one pass over the blocks with one task per
  (block, kind) that builds the index once, writes the zero-shift file,
  truncates and writes the one-shift file; the shortcuts are checked at run
  time against the managers' limits and fall back to the two stage
  processes. **The multiple-shift filter has no index family of its own**:
  its PRM diagonal index is the one-shift family's, the file
  `ptm_diag_index` (it used to be written twice, as `one_ptm_diag_index` and
  `multi_ptm_index`, byte-identical). `filter/mng/ptm_diag_index.hpp` holds
  that file name and the column limits (50000 / 100) from which both
  `OnePtmFilterMng` and `DiagFilterMng` initialize, so the two cannot
  diverge, and `OnePtmFilter` and `DiagFilter` check a file's column number
  against their limit after reading it (a mismatch means the index is built
  in memory with a warning). Keep it that way: a different limit for one of
  the two filters needs a file of its own. The joint tasks serialize
  **without** a mutex: each writes its own `MassMatch` into its own archive
  (Boost 1.90's singletons are built before `main`, so separate archives in
  separate threads do not interact), and a single writer would have been
  busy for most of the wall once every task writes two files.
  The filters are memory-bound: each thread streams every spectrum against
  one database block's index (117 blocks of ~1 GB one-shift index files for
  the human target+decoy database), so `compScores` misses cache on every
  increment and 16 threads scale about 2x, not 16x. In
  `massmatch/mass_match_util` the candidate finders take the index vectors
  by `const&` (they used to copy ~3.5 MB per call, four calls per spectrum
  set) and `findOneShiftTopProteinsWithRestrict` skips a protein whose best
  N-row plus best C-row score cannot reach the threshold and a row that
  cannot beat the current best; both are exact (the skipped pairs could not
  update the best) and cut the one-shift filter from 712 s to 389 s on
  480_F3_02. The scores are not sparse (about 40% of truncation rows score
  >= 1 per spectrum), so a sparse or sliding-window pair finder was slower
  and was not kept. The zero-shift filter has no pair scan (its finder is
  `simpleFindZeroShiftTopProteins` on per-proteoform maxima that
  `ZeroPtmFilter::computeBestMatch` computes once with
  `compProteoMaxScores` and shares between the four proteoform types); its
  time is the index scoring itself. The diagonal filter scores all peaks
  once per anchor peak (O(peaks^2 x index entries), ~40 ms per spectrum set
  on this database) with a full fill and max-scan of the rows per anchor;
  a touched-row list instead was 1.8x slower because most rows are touched.
- `src/search` — sequence/spectrum alignment (`diag`/`graph`/`graphalign`/
  `oneptmsearch`/`ptmsearch`/`varptmsearch`/`zeroptmsearch`). `graph` uses the
  header-only Boost Graph Library (`adjacency_list`, `graph_traits`, `graphviz`).
  All four search stages are multithreaded. The zero-shift and one-shift
  searches run one task per spectrum on a `SimpleThreadPool` through
  `common/thread/ordered_task_runner.hpp` (`OrderedTaskRunner<Result>`:
  results are handed back on the submitting thread in submission order, with
  a bounded number of tasks in flight), so they still write a single output
  file per proteoform type in spectrum order, byte-identical to the former
  sequential loops; the variable-PTM and PTM searches instead write one file
  per thread and merge them with `PrsmStrMerge`. Use the runner for a new
  stage whose output order matters; `ZeroPtmSearchMng` and `PtmSearchMng`
  carry the thread number from the toppic `-u` option.
- `src/stat` — E-value/p-value estimation (`count`/`local`/`mcmc`/`tdgf`),
  the PTM localization (`local`, and its alternative `local_align`) and
  the rescoring of the PrSMs before the FDR (`rescore`, `--rescore`);
  `tdgf/cluster_prot_corrector` is here rather than in `prsm` because it
  recomputes E-values (see `.claude/tools/toppic_output.md` for the
  localization and the correction, `.claude/tools/toppic_evalue.md` for
  the E-values and the rescoring).
- `src/merge` — the topdiff feature-merge backend (`feature_prsm`,
  `feature_sample_merge`), used by `console/topdiff_process`.
- `src/toplib` — toplib_build's spectral library building and
  toplib_search's library search (see `.claude/tools/toplib.md`); it
  depends only on
  `common` and `sql`.

These are all library layers. The **executable** layers are separate (they are
NOT in `COMMON_SRCS`; each is its own `add_executable` that links
`toppic_common`):

- `src/console` — the command-line tools. Each is a `main()` (`topfd.cpp`,
  `topdia.cpp`, `topindex.cpp`, `toppic.cpp`, `topmg.cpp`, `topdiff.cpp`,
  `toplib_build.cpp`, `toplib_search.cpp`; plus `topcalib.cpp`, a single
  file with its options in it, and the developer tool
  `theo_patt_gen.cpp`, a single file outside the default target) plus
  its `<tool>_argument.cpp` parser (`boost::program_options`) plus, for the four
  multi-step tools, its `<tool>_process.cpp` orchestrator. topfd/topdia have no
  `*_process` here — they drive the `topfd_process`/`topdia_process`
  orchestrators that live in the library. The `toppic_console_exe()` helper in
  `CMakeLists.txt` builds each, linking `toppic_common` + `Boost::program_options`.
- `src/gui` — Qt6 desktop front-ends (`topfd`/`topindex`/`toppic`/`topmg`/
  `topdiff`/`topdia`, plus `util` = a QProcess command builder + message
  helpers). The `toppic_gui_exe()` helper in `CMakeLists.txt` defines each
  target with per-target `AUTOMOC`/`AUTOUIC`/`AUTORCC` and
  `AUTOUIC_SEARCH_PATHS=src` (so the dialogs' `"gui/<tool>/ui_*.h"` includes
  resolve), linking `toppic_common` + `Qt6::Widgets/Core/Gui` +
  `Boost::program_options`. A dialog collects parameters and **shells out** to
  the matching CLI tool via QProcess, reading its default values from that
  tool's console argument parser. `src/gui/topmerge` is migrated but has no
  target (no `topmerge` CLI tool drives the `merge` backend yet).
  Every parameter control of the six forms has a tooltip (2026-10-12):
  rich text, `<html><head/><body><p>...</p></body></html>` (Qt wraps rich
  text; a plain-text tooltip is one long line), with `<`, `>` and `&` of
  the message HTML-escaped (`&lt;name&gt;`), the same text on a row's
  label and input, ending with the option it sets, e.g. `(topfd -c,
  --max-charge)`. Give a new control one, and change it with its option's
  help. The `-x` checkbox of toppic_gui and topmg_gui is "No TopFD feature
  file" ("Missing MS1 feature file" until 2026-10-12; it is the
  `<base>_ms2.feature` file, not an MS1 file). Start runs the form's `checkError()` first (every GUI since
  2026-10-12; before, only topmerge called it, so an empty field crashed
  topfd_gui and topdia_gui): the files (at least one, two for topdiff; no
  path with a space, `gui_check::pathWithSpace`, because the command line is
  split at spaces), and every line edit through `gui_check::invalidEdit`
  (empty, or not accepted by its validator), so a field's validator must
  accept exactly the values its option takes (the cutoff values take
  0.001 and 1e-5; thread numbers and ppm tolerances start at 1).

The `<tool>_argument.cpp` parsers are compiled into both the CLI tool and its
GUI (each target compiles its own object — no shared lib for them). When
migrating the remaining `*_process` drivers, drop the leftover
`xercesc::XMLPlatformUtils::Initialize()` calls — pugixml needs no global init.

When migrating a folder from the upstream Xerces tree, watch for include guards
that don't match the destination path (e.g. an `ms/env` file guarded
`TOPPIC_TOPFD_ENV_*`, or filter subdirs split as `ONE_PTM`/`VAR_PTM`) and
rename them to `TOPPIC_<PATH>_<FILE>_HPP_`. Three recurring upstream bites: a
missing `;` after a `LOG_ERROR(...)`/`LOG_DEBUG(...)` (our logger macro is a
`do {…} while(0)`, so it needs the terminator the upstream code omits); the
`str_util::toString(int)` calls that must become `std::to_string` (keep the
`toString(double)` ones — see the str_util table above); and a by-value Ptr/
container parameter that the function sorts/mutates in place, which must NOT be
"upgraded" to `const&` (it will not compile).
