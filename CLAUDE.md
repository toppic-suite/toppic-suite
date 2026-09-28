# Project conventions and migration notes

Guidance for working on this codebase (the `common/` and `seq/` layers of
TopPIC). Read this before adding or changing code in `src/common` or `src/seq`.

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
  the connection (opens it, creates the schema, closes it). When migrating the
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
- **ProteoWizard (pwiz)** — a trimmed copy is vendored under `ext/pwiz` (only the
  `utility/minimxml`, `utility/misc`, `data/common`, `data/msdata` source dirs
  are compiled; the rest is headers), built as a static `pwiz` library against
  Boost + zlib with `WITHOUT_MZ5` (so no HDF5). `ms/mzml/pw_ms_reader` uses it to
  read mzML/mzMLb. pwiz uses its own `minimxml` parser, **not** Xerces. It needs
  a few Boost pieces the system package lacks (`boost/nowide`, the `boost::enums`
  library, `foreach_field.hpp`); those are vendored under `ext/boost`, which
  `ext/` (the include root) resolves ahead of the system Boost while everything
  else still comes from the system. The vendored nowide's `iostream.cpp` (the
  Windows UTF-8 console streams `boost::nowide::cerr` etc. that pwiz's
  `Stream.hpp` uses) is appended to the `pwiz` sources; it is an empty TU off
  Windows, and it must be that copy, not the system `boost_nowide` library,
  because the vendored headers predate Boost.Nowide. Treat `ext/pwiz` and
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

- `install(TARGETS ${TOPPIC_EXECUTABLES} ...)` puts the six CLI tools and their
  six `*_gui` counterparts in `<prefix>/bin`. **When you add a new executable
  target, append it to the `TOPPIC_EXECUTABLES` list** (end of the file) so it
  is installed and gets the rpath below — nothing else is needed.
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
`toppic_manual.md`, `topmg_manual.md`, `topdiff_manual.md`. Each follows the
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
- TopDIA's `-p`/`-P` options set `TopfdPara::seed_env_inte_corr_tole_cutoff_`
  (copied into `EcscorePara` by `env_coll_detect::processMs1/processMs2`)
  before each feature-detection stage; both default to 0.5, the value the
  detector had hard-coded, so default runs are unchanged. TopPIC's hidden
  `-N`/`--n-terminal-label` works (its handler used to test the wrong key). (TopDiff resolves result
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

## Envelope refinement and the core-peak intensity fit (`ms/env/match_env_refine`)

After the DP stage, `DeconvSingleSp::postprocess` calls
`match_env_refine::mzRefine`, which regenerates each envelope's theoretical
distribution from its reference peak and rescales it to the experimental
peaks. It is two-stage, and the split is deliberate (issue 532):

1. **Monoisotopic choice** (current / previous / next isotope) uses the
   whole-envelope distance, unchanged from the original code. A core-only
   choice was tried and matched no more theoretical b/y fragments of horse
   myoglobin on a test MS/MS spectrum, so do not "simplify" stage 2's core
   into stage 1.
2. **Intensity ratio** of the chosen distribution is refit on the **core
   peaks only**: theoretical intensity >= `EnvPara::refine_core_ratio_`
   (0.5) x the reference peak, and present in the experimental envelope;
   fewer than 3 such peaks -> whole envelope (small envelopes unchanged).
   `compDistWithNorm`/`compDist` take an index list (empty = all peaks).
   Fitting every peak let tails inflated by overlapping neighbouring
   envelopes drag the ratio up, leaving theoretical apexes far above the
   observed ones (1.58x on the 16896 Da envelope near m/z 891 of
   `~/code/mms/data/myoglobin_peaks.txt`).

The theoretical intensities set here are what reach the outputs: the
`msalign` intensity per mass is the theoretical envelope's intensity sum,
the `.env` file's `THEO_INTE` column, and the SQLite `ms1_env`/`ms2_env`
rows. `assignIntensity` (called before `mzRefine` and again after sorting)
only splits an experimental peak among envelopes that share it; it never
changes theoretical intensities. When touching any of this, validate with
`topfd -T -g` on the myoglobin peak list: compare the deconvoluted masses
against the sequence's b/y ions (46 exact / 14 off-by-one-isotope is the
reference), and count envelopes whose theoretical apex exceeds the observed
apex by > 20% (77 of 611 before the core fit, 1 after).

## FLASHDeconv MS/MS deconvolution (`topfd/deconv/flashdeconv_ms2_process`)

`FlashDeconvMs2Process` is a second `DeconvMs2Process`: the same mzML group
reading, one task per MS/MS spectrum on the thread pool, per-thread msalign
files merged at the end and buffered SQLite records, but the spectrum is
deconvoluted with the FLASHDeconv algorithm (OpenMS `SpectralDeconvolution`,
`PeakGroup`, `PeakGroupScoring`, `FLASHHelperClasses`) instead of
`DeconvSingleSp`. The option `topfd --ms2-flashdeconv` selects it
(`TopfdPara::isUseFlashDeconv`, `topfd_process.cpp`). Facts to keep straight:

- **Everything lives in the cpp**, in the namespace `flashdeconv_ms2_process`:
  `Averagine` (the theoretical distributions), `LogMzPeak`, `PeakGroup` (one
  deconvoluted mass with its scores) and `FlashDeconvSingleSp` (one spectrum).
  Each function's comment names the OpenMS function it ports; keep that when
  editing so the two can be compared.
- **FLASHDeconv's constants, not TopFD's.** The tolerance is 10 ppm
  (`kMs2TolerancePpm`), the isotope spacing 1.002371 Da (`kIsoDaDistance`, not
  1.00335), the log m/z bin is tolerance / 2.5, the thresholds are cosine 0.85,
  charge SNR 0.25, QScore 0.2, charge score 0.7, and the QScore weights are
  the fixed ones of `PeakGroupScoring`. TopFD's `-e` (m/z error) is used only
  to match the theoretical envelopes to the peaks for the output.
- **The averagine table is EnvBase's** (`res/base_data/theo_patt.txt`, one
  entry per 25 Da from 50 Da to `-m`), prepared as `PrecalculatedAveragine`
  prepares CoarseIsotopePatternGenerator's: the weak end isotopes are trimmed
  until 0.01% of the power is gone, the intensities are L2-normalised (so
  `getCosine` needs no norm of the theoretical side), and the left / right
  isotope counts from the apex are a running maximum over the masses. It is
  built once per run in `process()` and shared read-only by the threads.
- **Output.** Every surviving peak group becomes a `MatchEnv` of its
  representative (highest per-charge SNR) charge: the EnvBase envelope of the
  monoisotopic mass, scaled so its intensity sum is the group's intensity over
  all charges, **trimmed to the noise floor exactly as `env_detect::compFullEnv`
  trims the standard pipeline's own envelope** (`getSubEnv(1.0, bound_inte,
  INT_MAX, INT_MAX)`; `bound_inte`/`min_inte` are `deconvMsTwo`'s own estimate
  of `DeconvData::initMinInte`, via `baseline_util::getBaseLine` on the
  spectrum, gated by `-estimate-min-inte`/`isEstimateMinInte()`) before
  matching to the peaks by `ExpEnv`. Without this trim every envelope carried
  the averagine table's full isotope tail (EnvBase's entries run to 100+
  peaks for large masses) into the SQLite output, 3–4x more than
  `DeconvMs2Process` writes for the same data (toppic_test experiment,
  2026-09). The QScore is stored as the EnvCNN score (what
  `match_env_util::getDeconvMsPtr` reports and sorts by) and the isotope
  cosine as the MsDeconv score; no EnvCNN model runs. Note FLASHDeconv's own
  msalign writer reports the maximum charge of the range, this path the
  representative charge.
- **Not ported:** the MS1 side and the precursor peak group search (the MS/MS
  max mass and charge come from TopFD's precursor features, as in
  `DeconvMs2Process`), decoy / FDR runs, target masses, and the tolerance
  auto-estimation.

## Post mass matching (`prsm/prsm_post_mass_match`, toppic only)

After the search, the clustering, the cluster protein correction and the
match recount, `console/toppic_process` runs
`prsm_post_mass_match::process`, which matches the theoretical fragment
masses of each PrSM that no deconvoluted mass matched against the
**centroided MS/MS peaks stored in TopFD's SQLite database** (MSPathFinderT's
envelope test), scores the added masses with EnvCNN, and writes
`<base>_post_ms2.msalign` plus PrSM/feature files under that name, which every
downstream step (PTM characterization, FDR, tables, SQL output) then reads
instead of the input. Facts to keep straight:

- **It is on by default.** `initArguments()` sets `postMassMatch = "true"`;
  `-D` / `--disable-post-match` sets it to false (there is no enabling flag
  any more), and `-I` / `--post-min-peak-num` (default 1) is the minimum
  number of observed isotopic peaks per accepted fragment.
- **A missing database is an error, not a skip.** The database name comes
  from `prsm_post_mass_match::sqlFileName(sp_file)` (`<base>.sqlite` for
  `<base>_ms2.msalign`; use it, do not re-derive the name). It is checked in
  three places with the same message ("run TopFD without disabling its SQLite
  database output (without -N) ..., or run TopPIC with --disable-post-match"):
  `ToppicArgument::validateArguments` (so a run fails before the search),
  `process()` itself (abort), and `toppic_gui`'s `checkError()` (a warning
  dialog before toppic is launched). Keep all three in step.
- **It runs on the best PrSM of every spectrum, after the clustering, the
  cluster protein correction and the recount.** `TopPIC_postprocess`
  clusters and corrects the PrSMs (`toppic_raw_prsm` → `toppic_cluster_prot`),
  recounts them against the full spectra (`toppic_recount`), then
  post-matches `toppic_recount` (targets and decoys alike) into
  `<base>_post_ms2.toppic_recount`,
  and the PTM characterization (`-l`), FDR, cutoffs, tables and SQL output
  read the spectra with the added masses under the new name. Post matching
  changes neither E-values nor cluster ids, so the FDRs and the
  identifications equal those of a `-D` run; only the matched masses and
  fragments differ. (It used to run before the clustering, which let the
  cluster protein correction recompute E-values on the post-matched
  spectra: 42 E-values and 6 PrSMs at the 1% cutoff differed from `-D` on
  480_F3_02, toppic_test report 0062.) The PrSM XML carries the cluster
  ids (`proteo_cluster_id`, `prot_id`, `prot_cluster_id` in the proteoform
  element), so the read-modify-write in `process()` keeps them.
- **It reads each PrSM's peaks with one query on `ms2_peak` by `spec_id`**
  (`RawPeakReader`), and TopFD writes that table without an index, so the
  reader first runs `CREATE INDEX IF NOT EXISTS ms2_peak_spec_id_index`.
  Without the index every query scans the whole table: on 480_F3_02
  (4.5 M peaks, ~4,000 PrSMs) the step took 11 of the 18 minutes of a
  default run. Keep any new per-spectrum query on the SQLite file indexed.
- **Cluster protein correction runs right after the clustering, before
  the match recount** (`stat/tdgf/cluster_prot_corrector`; it lives in
  `stat` because it needs `EValueProcessor`): in every proteoform cluster
  the PrSMs whose protein differs from the protein of the cluster's
  lowest-E-value PrSM (the same proteoform mass matched to another
  protein) get that PrSM's proteoform (shared pointer; its mass becomes
  their adjusted precursor mass), and their matched fragment numbers and
  E-values (`EValueProcessor::compEvalues`) are recomputed on the spectra
  read as the search read them: `MsAlignReader` with the EnvCNN cutoff and
  `SpPara::getFragDedupPpo()` (the combined masses of `--frag-mass-dedup`);
  the recount that follows replaces the matched numbers with full-spectrum
  counts (a `TdgfMng` built from the same values as the search's, with
  `varPtmTypeNum` from the argument map, 0 if absent, and, with
  `--topmg-identify`, max(`-s`, `--topmg-num-shift`) unexpected shifts, since
  a graph-alignment proteoform is corrected like any other). A corrected PrSM
  also takes the best PrSM's `<tool>` (toppic or topmg), so the per-search
  FDR counts it with that search. Only with `--topmg-use-mcmc` are the other
  proteins' PrSMs of a cluster whose best PrSM comes from the graph alignment
  dropped instead (their generating-function E-value would be counted
  against MCMC ones; `Cluster protein correction - removed N PrSMs`).
  Otherwise nothing is removed; a corrected PrSM that matches at most
  `comp_evalue_min_match_theo_frag_num_` fragments gets the maximum
  E-value and fails the cutoffs. The output is `toppic_cluster_prot`
  (deleted without `-k`); the recount reads it and writes `toppic_recount`
  with the matched numbers of every PrSM, corrected ones included, counted
  on the full spectra, which the post mass matching, the PTM
  characterization (`-l`) and the FDR use.
  Combined runs (`TopPIC_multi_postprocess`) do not run it. The earlier
  `prsm_cluster_prot_filter`, which dropped those PrSMs, is gone.
- **`<base>.toppic_record_prsm` is where `TopPIC_postprocess` ends and
  `TopPIC_report` starts.** It is a copy, made inline at the end of
  `TopPIC_postprocess`, of the PrSM file after the FDR computation
  (`toppic_cluster_fdr` with `-d`; `toppic_recount`, or
  `toppic_cluster_local` with `-l`, without): the best PrSM of every
  spectrum with the added masses, the cluster ids and the FDRs, before the
  cutoffs, **decoys included** (`prsm_fdr::process` and
  `prsm_fdr_groups::process` always write the decoys; their former
  keep-decoy parameter is gone). `cleanToppicDir` never deletes it. `TopPIC_report`
  (the decoy removal, the cutoffs, the tables and the SQL output) reads
  only the argument map and re-derives the file names from it, so the two
  stages can be run separately, and a combined run (`-c`) merges the
  fractions' record files and runs the report alone. Without `-K`,
  `prsm_decoy_filter::process` writes the PrSMs whose sequence name does
  not start with `DECOY_` to `toppic_prsm_target`, which the spectrum-level
  cutoff reads (with `-K` it reads `toppic_record_prsm` itself); the
  proteoform and protein cutoffs and their best-PrSM selections chain from
  `toppic_prsm_cutoff`, and `cleanToppicDir` deletes `toppic_prsm_target`
  and the cutoff files with the other intermediates.
- **The PrSM-only output steps take the spectrum base name, not the
  spectrum file.** `prsm_decoy_filter`, `prsm_cutoff_selector`,
  `prsm_form_filter`, `prsm_prot_filter` and `prsm_sql_output::write` only
  ever name `<base>.<ext>` files (the SQL writer strips `_ms2`/`_post_ms2`
  from it for the `.sqlite` name), so their first name argument is
  `base_name`, the spectrum file name without its extension
  (`file_util::basename(sp_file)`); `TopPIC_report` and `TopMG_post`
  compute it once. Do not pass the spectrum file name to them. The steps
  that need the spectra themselves — `PrsmMatchTableWriter` through
  `prsm_reader_util::readPrsmsWithSpectra`, and `PrsmCoverage` — open the
  msalign named in `PrsmPara`, so the report still needs `<base>.msalign`
  on disk for the tables (the PrSM XML stores no peaks).
- **Result file names carry `_post_ms2`** whenever post matching ran; the
  cleanup in `toppic_process` (`cleanToppicDir`) takes that name to find the
  intermediate files. With `-D` the names are the plain `<base>_ms2` ones.
- **Combined runs (`-c`) merge the spectra and record files, re-cluster,
  compute one FDR and report.** `TopPICProgress_multi_file` merges each
  fraction's post-matched spectra (`<frac>_post_ms2.msalign`, or
  `<frac>_ms2.msalign` with `-D`) into `<name>_ms2.msalign` with
  `msalign_frac_merge::mergeFractions`, and each fraction's
  `toppic_record_prsm` (`<frac>_post_ms2.toppic_record_prsm`, or
  `<frac>_ms2.…` with `-D`) into `<name>_ms2.toppic_raw_prsm`; `<frac>` is
  the mzML base without the `_mono` suffix of a `-g` run (the same rule as
  `prsm_post_mass_match::sqlFileName`). The spectra are merged from the
  post-matched files because the report attaches them to the PrSMs, whose
  matches were counted on those spectra. `TopPIC_multi_postprocess` then
  clusters the merged PrSMs into proteoforms and proteins across the
  fractions with `prsm_simple_cluster` (by protein and precursor mass, not
  by feature; the feature files are not merged) and computes one FDR
  (decoys kept), writing
  `<name>_ms2.toppic_record_prsm`, and `TopPIC_report` runs on it. The
  re-clustering is required: `prsm_util::mergePrsmFiles` offsets the
  spectrum and feature ids per fraction but not the proteoform/protein
  cluster ids, so without it the fractions' clusters collide and the
  best-PrSM filters collapse them (79 combined proteoforms instead of 142
  on st_1 + st_2). There is no recount, post mass matching or PTM
  characterization on the merged set. The record file is always kept, so
  `-C` (combine-result-only) does not need `-k`. Combined results are
  always `<name>_ms2_toppic_*`.
- **The multiple-shift stage runs after the first-pass top selection
  (`TopPIC_identify` → `TopPIC_multi_ptm_search` → `TopPIC_postprocess` →
  `TopPIC_report`, no `TopPIC_post` any more; a combined run's merged set
  goes `TopPIC_multi_postprocess` → `TopPIC_report`).** `TopPIC_multi_ptm_search`
  reads `toppic_pre_raw_prsm` (the first-pass top selection), skips the
  spectra whose best PrSM has an E-value strictly below the skip E-value
  (with `-d` the smallest E-value of a decoy PrSM of the first pass, so
  that every spectrum a decoy could compete with is searched again; 0.01
  without a decoy database; the hidden `--multi-skip-evalue <x>` overrides
  both), searches the rest (`DiagFilterMng::skip_spec_ids_`), merges all
  first-pass PrSMs (`toppic_pre_raw_prsm`, skipped spectra included) with
  the multiple-shift ones and keeps the best PrSM per spectrum
  (`toppic_raw_prsm`; without `-s 2` `TopPIC_postprocess` copies
  `toppic_pre_raw_prsm` to that name). `TopPIC_postprocess` runs the
  clustering, cluster protein correction, recount, post mass matching,
  characterization and FDR once on `toppic_raw_prsm` and writes
  `toppic_record_prsm`;
  `TopPIC_report` runs the decoy removal, the cutoffs and the outputs.
  The PrSM file chain is therefore `toppic_pre_raw_prsm` → `toppic_raw_prsm`
  → `toppic_proteoform_cluster` → `toppic_cluster` → `toppic_cluster_prot`
  → `toppic_recount` [→ the same under the `_post_ms2` name]
  [→ `toppic_cluster_local`] [→ `toppic_cluster_fdr`]
  → `toppic_record_prsm`
  [→ `toppic_prsm_target`] → `toppic_prsm_cutoff` → `toppic_form_cutoff`
  → `toppic_form_cutoff_form` / `toppic_prot_cutoff` → `toppic_prot_cutoff_prot`.
  Do not reintroduce a separate FDR for the multiple-shift PrSMs: with the
  skipped targets missing from its denominator it under-reports (toppic_test
  reports 0023/0024).
- **GUI wiring.** `toppicwindow.ui` has `postMassMatchCheckBox` (checked by
  default from the parser's argument map, enables `postMinPeakNumEdit`) on the
  bottom row of the Advanced Parameters tab; `gui/util/command.cpp` emits
  `-D` when the box is **cleared** and `-I <n>` only while it is checked.
- **Short options in `toppic_argument.cpp` are nearly exhausted.** Before
  adding one, list the letters in use with
  `grep -oE '"[a-z-]+,[A-Za-z]"' src/console/toppic_argument.cpp`; a clash is
  only reported at run time as "option '-X' is ambiguous" (that is how `-P`
  collided with `--proteoform-ppm-error`). Every option is declared twice
  (the visible `display_desc` and the full parse `desc`); change both.

## TopDiff (`merge/feature_sample_merge`, `console/topdiff_process`)

TopDiff reads one `<base>_post_ms2_<tool>_proteoform.xml` (or `_ms2_`)
per sample into `FeaturePrsm`s (protein, residue range, adjusted precursor
mass, proteoform cluster id and abundance, feature min/max/apex retention
times), normalizes each sample's times by its latest feature end, aligns
samples 2..n onto sample 1 (`sampleAlignTime`: a dynamic programme over
the 1,000 most abundant proteoforms of each, sorted by apex time; a cell
matches when a feature of the other sample within 0.1 normalized time has
a precursor mass within `--error-tolerance`, scored by log10 abundance;
the backtracked matches are the piecewise-linear warp `setAlignTime`
applies), then builds the table (`getFeatureTable`): all proteoforms of
all samples by decreasing abundance, each not yet in the table starts a
row and takes, in every other sample, the most abundant unused proteoform
of the same protein within `--error-tolerance` and 0.3 aligned apex time
(`findMatchFeature`). Facts to keep straight:

- **The matching is by protein accession, precursor mass and aligned
  time only**; the residue range and the proteoform string are not
  compared. On 480_F3_02 result sets of the same TopFD features
  (toppic_test experiment 65) that matches every proteoform of a sample
  with itself and the same spectrum (1,765 of 1,765), and across search
  settings 1,780 of 1,799 matched rows join the same MS/MS scan; 4 rows
  join different features of one protein within 1.2 Da.
- **The hot spots were exact rewrites** (experiment 65, byte-identical
  tables): the alignment's per-cell `findMatch` scan became per-row sorted
  match-position lists with a binary search (`matchPositions` /
  `findMatch`; three quarters of the run before), the matrix is a flat
  `std::vector<Cell>` with the cell matches kept for the backtracking,
  `featureUsed` is a set of (sample, cluster id), and `findMatchFeature`
  looks only at the protein's features through a per-sample index in the
  intensity order of the scan it replaced (the dead `prot == ""` branch is
  kept, with the strict tolerances, so the result is the same). 32 copies
  of one result set: 14.8 s -> 2.4 s; two samples: 0.56 s -> 0.17 s, most
  of it now the XML reading. `initCells` used to leave cell (0, b_size)
  uninitialized (fixed; no output change on that data).
- **Result files of other versions are refused**, not parsed: every
  `<prsm>` starts with a `<version>` element (`Prsm::toXmlElement`,
  `Version::getVersion()`), and `PrsmReader::readOnePrsmStr` exits with
  "was written by version X of TopPIC Suite and cannot be read by this
  version" when it differs from the running version, or "by an earlier
  version" when it is absent (1.9 and earlier wrote `match_peak_num` /
  `match_fragment_num` and no `rt`; `PrsmStr` threw `std::invalid_argument`
  on them). Every tool that reads PrSM files gets this check; a result
  file needs the search rerun after a version change.
- The console prints `Aligning retention times of sample <k>` and
  `Matching proteoforms across samples` as started/finished steps (the
  manual lists them); the table is written next to the **first** input
  file, so give inputs from one directory or expect the output there.
- Not done: matching an identified proteoform to an unidentified TopFD
  feature of another sample (the commented-out "Match without MS/MS ID"
  path) - the `_ms1.feature` files are not read.

## Average masses in msalign files (`ms/env/avg_mass_util`)

`topfd --average-mass` (`-A`, `TopfdPara::output_avg_mass_`) writes the
MS/MS msalign file with average masses: `deconv_ms2_process::deconvMsTwo`
and `topfd_single_process` pass the `DeconvMs` through
`avg_mass_util::getAvgMassMs` just before `MsAlignWriter::writeMs`, which
copies the header (precursor m/z and mass) and the peaks with every
monoisotopic mass replaced by `EnvBase::convertMonoMassToAvgMass` (the
intensity-weighted average of the `theo_patt.txt` distribution, cached in
`Env::avg_mz_`). The SQLite record, the feature files and the MS1 file keep
monoisotopic masses. `toppic --average-mass` (`-g`, argument key
`averageMass`) converts each input file in `TopPICProgress_multi_file`,
before any search, with `avg_mass_util::convertToMonoMass` to
`<base>_mono_ms2.msalign` (`EnvBase::convertAvgMassToMonoMass` inverts the
table by fixed-point iteration; the feature file is copied to the new
name) and replaces the spectrum file list with the converted names, so the
per-file runs, the cleanup and the `-c` merge all see monoisotopic files.
`prsm_post_mass_match::sqlFileName` strips the `_mono` suffix so the post
mass matching finds `<base>.sqlite` and still writes `<base>_post_ms2.*`.
Do not convert masses inside `MsAlignReader` instead: the converted file
is what every later step (and the user) sees. GUI: `toppicwindow.ui` has
`averageMassCheckBox` (Advanced Parameters tab, above the localization
box; the keep-decoy row moved up to make room), read into
`arguments_["averageMass"]`, and `gui/util/command.cpp` emits `-g` when it
is checked.

## Precursor and fragment mass tolerances (`para/sp_para`)

`SpPara` holds two `PeakTolerance` objects, `prec_tolerance_ptr_` and
`frag_tolerance_ptr_`, set by `-E`/`--prec-error-tolerance` and
`-e`/`--frag-error-tolerance` (argument keys `precErrorTolerance` and
`fragErrorTolerance`, both default 10 ppm) in toppic, topindex and topmg and
their GUIs. The index file name (`filter/mng/index_file_name`) carries only
the fragment value: the indexes depend on neither tolerance, and the manuals
say an index is reused when `-f`, `-n`, `-e` and `-d` match. topmg's hidden
`--proteoform-ppm-error` moved from `-E` to `-Q` to free the letter. The
rule for choosing between the two tolerances: the
precursor tolerance covers anything sized on the precursor mass — the
`getFirstPrecErrorTolerance(ppo)` checks in the filters/searches/E-value
code, `refinePrecursorAndHeaderShift`, the localization shift checks
(`LocalMng::prec_tole_ptr_`), **and** the tolerance handed to
`prm_ms_util::getIntMassErrorList` / `getPrmPeakPtrs` and
`comp_prob_value`, which use it only to build the synthetic zero and
precursor PRM peaks (`prm_peak_factory`, sized on the precursor mass; the
real fragment peaks already carry their own tolerances). The fragment
tolerance covers per-peak tolerances set on fragment masses: the PRM /
extend spectrum factories, post mass matching's m/z test, and the
localization M tables (`LocalMng::frag_tole_ptr_`). Neither `SpPara` nor
`PeakTolerance` has XML (de)serialization any more: their `appendXml` and
`XmlDOMElement` constructors were never called and were removed; `SpPara`
is built only from the toppic/topmg arguments in `PrsmPara`.

## Precursor +-1 Da error per identification stage (`para/sp_para`)

`SpPara` holds one precursor-offset list per stage; each filter/search
expands a spectrum into one spectrum set per offset with
`spectrum_set_factory::geneSpectrumSetPtrVecWithPrecError`. Defaults: zero
shift and one shift `{0, -1 Da, +1 Da}`, variable PTM `{0}`, multi shift
`{0}` (the diagonal filter and PTM search absorb a 1 Da error into the
unknown shift, so they have no switch). `PrsmPara` applies the toppic hidden
options `--disable-zero-shift-prec-error`, `--enable-one-shift-prec-error`
and `--enable-var-ptm-prec-error` (argument keys `zeroShiftPrecOneDaError`,
`oneShiftPrecOneDaError`, `varPtmPrecOneDaError`) through the
`SpPara::set*PrecOneDaError` setters. **toppic's defaults differ from
`SpPara`'s for one stage**: `ToppicArgument::initArguments` sets
`oneShiftPrecOneDaError = "false"` (one-shift stage searched with `{0}`;
`--enable-one-shift-prec-error` restores `{0, -1, +1}`) and
`varPtmPrecOneDaError = "false"` (variable-PTM stage searched with `{0}`,
matching the `SpPara` default; `--enable-var-ptm-prec-error` switches it to
`{0, -1, +1}`, the pre-2026-09 toppic default). Absent keys
keep the `SpPara` defaults, so topmg is unaffected. The fragment-mass
`frag_ext_offsets_` expansion (`{0, -1 Da, +1 Da}` applied by
`extend_ms_factory` to every fragment above `frag_extend_min_mass_`) is separate:
`SpPara::setFragMassExtension(false)` reduces it to `{0}`, reached through
toppic's hidden `--disable-frag-mass-extension` (key `fragMassExtension`,
default `"true"`), and `setFragExtendMinMass` sets the threshold, reached
through the hidden `--frag-extend-min-mass <Da>` (key `fragExtendMinMass`,
default `"5000"`, the `SpPara` default; validated as a non-negative
number).

## Source layout

`src` is the include root, so headers are included by their path from `src`
(`common/...`, `seq/...`, `ms/spec/...`, `para/...`, `sql/...`). All of the
following compile into the single `toppic_common` **OBJECT** library (the
`COMMON_SRCS` glob in `CMakeLists.txt` lists each directory) — its object files
are linked directly into each executable, so **no `libtoppic_common.so`/`.a` is
built**; the executables are self-contained. Each layer depends only on the ones
above it, never the reverse:

- `src/common` — foundation: `base`, `util`, `xml`, `thread`.
- `src/sql` — thin SQLite helper (`sql_util`: `execSql`, `prepareSql`,
  `stepAndReset`, all abort on error — use them rather than re-implementing
  the prepare/step boilerplate) plus `ms_sql_writer`, the
  MS1-peak tables for 3D visualization: `MsSqlWriter` writes a peak list as `PEAKS0`
  (all peaks) plus down-sampled `PEAKS1..n` layer tables, one `CONFIG` row
  per layer, and a colour bucket per peak, on a caller-owned connection like
  the topfd writer below. topfd `--sql-3d` (in `deconv_ms1_process`) writes
  them into the topfd `.sqlite` after MS1 deconvolution.
- `src/para` — analysis parameters (`sp_para`, `peak_tolerance`, `prsm_para`).
- `src/seq` — sequence/proteoform layer.
- `src/ms` — mass-spectrum layer: `spec` (peaks/spectra/msalign), `msmap`,
  `factory`, `env` (envelope detection), `feature`, `mzml`.
  `MsAlignReader` has a four-argument constructor taking an EnvCNN score
  cutoff (`SpPara::getEnvCnnCutoff()`, whose own default is 0.2, but both
  toppic and topmg set the key from their `-F`, which defaults to 0.4 in
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
  line, space-separated) with `DeconvSingleSp`, writing `<input>_ms2.msalign`,
  `<input>_ms2.env` and `<input>.sqlite`. **Both** name their outputs through
  `TopfdPara::setMzmlFileNameAndFaims` (base name = input path minus
  extension, `_<voltage>` appended for FAIMS; it also creates the SQLite
  database when enabled) — do not hand-roll output names or a second
  `createSqlDb` call in a driver.
- `src/topdia` — the TopDIA pseudo-spectrum layer built on `topfd`/`ms`:
  `common` (`topdia_para` + `topdia_process`) and `pseudo_spec`
  (`mzrt_feature`, `pseudo_peak`, `pseudo_spectrum`, `generate_pseudo_spectrum`).
- `src/prsm` — proteoform-spectrum-match layer (`Prsm`/`SimplePrsm`, readers,
  writers, FDR, clustering, coverage). `prsm_reader_util::readPrsmsWithSpectra`
  reads a PrSM file and attaches the deconvoluted/refined spectra from the
  msalign file; both `prsm_match_table_writer` (the TSV tables) and
  `prsm_sql_writer` (the same columns as `prsm`/`proteoform`/`protein` tables, plus
  `prsm_mass_shift` and `prsm_protein_match`) are built on it.
  `prsm_sql_output::write` is the one entry point that writes those tables
  plus `seq/fasta_sql_writer`'s `fasta_seq` table into the topfd `.sqlite` of
  the spectrum file when it exists; toppic and topmg both call it after their
  protein table, passing their own cutoff file extensions — do not duplicate
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
- `src/stat` — E-value/p-value estimation (`count`/`local`/`mcmc`/`tdgf`);
  `tdgf/cluster_prot_corrector` (see the post mass matching section) is
  here rather than in `prsm` because it recomputes E-values.
- `src/merge` — the topdiff feature-merge backend (`feature_prsm`,
  `feature_sample_merge`), used by `console/topdiff_process`.

These are all library layers. The **executable** layers are separate (they are
NOT in `COMMON_SRCS`; each is its own `add_executable` that links
`toppic_common`):

- `src/console` — the command-line tools. Each is a `main()` (`topfd.cpp`,
  `topdia.cpp`, `topindex.cpp`, `toppic.cpp`, `topmg.cpp`, `topdiff.cpp`) plus
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
