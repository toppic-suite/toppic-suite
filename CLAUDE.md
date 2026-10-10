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
- TopDIA's `-p`/`-P` options set `TopfdPara::seed_env_inte_corr_tole_cutoff_`
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

## Quantized ONNX models (`res/envcnn_models`, `res/ecscore_models`)

The models the tools load, `envcnn_two_block.onnx` and
`ecscore_seven_attr.onnx` (file names fixed in `onnx_env_cnn::initModel` and
`onnx_ecscore::initModel`), are dynamically quantized (int8 weights,
activations quantized at run time) forms of the float PyTorch exports kept
next to them as `envcnn_two_block_full.onnx` and
`ecscore_seven_attr_full.onnx`, which no code loads (toppic_test experiment
131, issue #595). Both are made by `quantize.py` of that experiment (Python
onnxruntime 1.24.1; the files keep opset 11 and IR version 6 and load in the
vendored runtime 1.23.2). They are used by TopFD and TopDIA (both models) and
by TopPIC's post mass matching (EnvCNN).

- **Gemm -> MatMul + Add first.** The export writes every fully connected
  layer as `Gemm(A, W, b, transB = 1)`, which onnxruntime's
  `quantize_dynamic` does not quantize; the rewrite `MatMul(A, W^T) + b`
  gives the same logits.
- **EnvCNN (`fc_row`).** The three fully connected layers (fc.0 9,600 ->
  2,048 holds 19.7 M of the 21.9 M parameters) get int8 weights with one
  scale per output channel, and their inputs are quantized to uint8 with
  **one scale per row** (envelope): `max(row) / 255`, valid because each is
  a ReLU output (>= 0), built from ReduceMax, Max, Div, Round, Cast,
  MatMulInteger and Mul. `quantize_dynamic` itself uses one activation scale
  per call: the deconvolution scores up to 256 envelopes per call, and its
  input has a few rows with values up to 825 (experimental intensities
  divided by the theoretical envelope's highest peak; a neighbouring intense
  peak), which coarsened the whole batch: scores moved by 0.06 on average,
  4.7% of the fragment masses crossed toppic's 0.4 (`-F`), and 142 PrSMs of
  480_F3_02 changed their proteoform. With one scale per row a score does not
  depend on the batch (mean change 0.002). The convolutions stay float:
  ConvInteger was slower than the float convolutions and less accurate. The
  file is 22 MB instead of 87 MB.
- **ECScore (`fc`).** `quantize_dynamic` with int8 weights, one scale per
  tensor (113 k parameters, 122 kB). TopFD calls it once per candidate
  feature (about 1 s per run in total either way); its scores change by
  0.027 on average (its seven inputs have different ranges but share one
  activation scale).
- **Results** on 480_F3_02 (toppic_test experiment 131, report 0131; main
  7eb5dbe9 with the models swapped in `res/`, i7-11700 with AVX-512 VNNI,
  `topfd -u 16`, three runs each): TopFD 283 s against 408 s with the float
  models (MS1 deconvolution 144 against 216 s, feature detection 64 against
  75 s, MS/MS deconvolution 74 against 115 s; EnvCNN was 84-88% of the
  deconvolution time and its int8 calls are 1.6-1.8x faster). 42,447
  features against 45,360 (the quantized ECScore; the EnvCNN alone gives
  45,367), 320 of the 4,318 MS/MS precursors changed. toppic (setting d of
  toppic_test report 0062): 3,346-3,347 PrSMs, 1,996-1,998 proteoforms and
  436 proteins against 3,336-3,342 / 1,998-2,000 / 437; 33 PrSMs with
  another proteoform than in a float run (1-3 between float runs). The int8
  speed-up depends on the CPU's int8 dot-product instructions; it was not
  timed on CPUs without AVX-512 VNNI or on ARM.
- **In the tree** (branch `quantized_onnx_models` on main 72534ded, merged
  into main on 2026-10-06; toppic_test experiment 131, amendment 1; one run
  each): TopFD 278.6 s
  against 402.8 s with the float `*_full` models; its outputs are as close to
  the experiment's arm B as two runs of B are to each other. toppic: 3,347
  PrSMs, 1,997 proteoforms, 437 proteins against main's 3,338 / 2,001 / 436.
  TopPIC's post mass matching with the int8 EnvCNN (5.7 s of EnvCNN against
  16.9 s) gives the same tables and XML files as with the float model on the
  same TopFD output; only 3,772 of the 724,767 EnvCNN scores in
  `_post_ms2.msalign` change (up to 0.12, 41 across 0.4). TopDIA was not
  tested.
- **The mass-difference input of unmatched peaks** (`onnx_env_cnn::
  extractFeature`, and ECScore's EnvCNN feature in `comp_env_cnn_score`): a
  theoretical peak without an experimental peak within 0.02 m/z gets -0.02
  or +0.02, alternating over the peaks (even index -0.02). It was a random
  sign from `std::rand() % 2`, one stream for the process, so an envelope's
  score depended on how many envelopes had been built before it and, with
  several threads, on their scheduling: the "EnvCNN jitter" of two runs of
  one binary (toppic_test report 0036 named the cause). toppic_test
  experiment 139: a PTCR spectrum gave 95 or 94 species at `topfd -u 1`
  depending on whether an HCD scan was deconvoluted before it in the same
  thread (scores 0.562 against 0.561 moved an envelope across 0.2). With
  the alternating sign two `topfd -u 16` runs on 480_F3_02 write the same
  `_ms1.msalign`, `_ms2.msalign` and feature files byte for byte and the
  same SQLite rows (in another insertion order; each thread buffers its
  records); two runs of main 6eb40e97 differ in 771 of 2,803 MS1 spectra
  and 193 MS1 features. Against main the output lies within main's own
  variation (814 MS1 spectra, 131 MS1 features), and toppic in setting d
  of toppic_test report 0062 gives 3,355 PrSMs, 1,998 proteoforms, 438
  proteins, main's first run's 3,355 scans and proteins (its second run
  3,354 / 1,997 / 438). TopPIC's post mass matching uses the same input
  builder, so its EnvCNN scores are repeatable now too (not measured).

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
changes theoretical intensities. The greedy MS1 and MS/MS searches (TopFD's
defaults) do not run this refinement. When touching any of this, validate
with `topfd -T -g --deconv-ms2 msdeconv` on the myoglobin peak list (without
the option `-T` runs the greedy search): compare the deconvoluted masses
against the sequence's b/y ions (46 exact / 14 off-by-one-isotope is the
reference), and count envelopes whose theoretical apex exceeds the observed
apex by > 20% (77 of 611 before the core fit, 1 after).

## Deconvolution methods (`--deconv-ms1`, `--deconv-ms2`)

topfd and topdia choose the deconvolution of each spectrum level with two
visible options: `--deconv-ms1 <greedy|msdeconv>` and `--deconv-ms2
<greedy|msdeconv|flashdeconv>` (`greedy`: the greedy searches of the next
sections; `msdeconv`: MS-Deconv's dynamic programming, `DeconvSingleSp`;
`flashdeconv`: the FLASHDeconv port below). Both tools also take the
short names `-p` and `-q` (2026-10-08; in topfd briefly `-x` and `-y`); in
topdia `--ms1-intensity-correlation-cutoff` moved from `-p` to `-y` for
them and `--ms2-intensity-correlation-cutoff` from `-P` to `-Y` (topdia_gui's
command emits `-y` and `-Y`). Facts to keep straight:

- **`TopfdPara` stores the names** (`ms1_deconv_method_`,
  `ms2_deconv_method_`; `getMs1DeconvMethods()` / `getMs2DeconvMethods()`
  list the accepted values, `setMs1DeconvMethod` / `setMs2DeconvMethod`
  return false for any other, which both parsers report as an error). The
  drivers ask `isMs1GreedyDeconv()`, `isMs2GreedyDeconv()` and
  `isUseFlashDeconv()`, so a new method is a new name plus its branch in
  `deconv_ms1_process::geneTask`, `deconv_ms2_process::deconvMsTwo`,
  `topfd_single_process` or the `topfd_process` / `topdia_process` choice
  of the MS/MS process class.
- **Defaults differ by tool**: `TopfdPara`'s own default, TopFD's, is
  `greedy` for both; TopDIA's `getTopfdParaPtrForTopdia` sets `msdeconv` for
  both (the greedy searches were measured on DDA data only), so a default
  TopDIA run is unchanged. TopDIA's MS/MS step (`ms2_raw.msalign`, run with
  `setMissingLevelOne(true)`) can use `FlashDeconvMs2Process`, which handles
  the missing-MS1 mode as `DeconvMs2Process` does.
- **`-T` refuses `flashdeconv`** (the FLASHDeconv process reads mzML files
  only); before, `-F` was silently ignored there.
- **Parameter block**: "MS1 deconvolution: Greedy" and "MS/MS
  deconvolution: Greedy" / "FLASHDeconv" are printed only for those methods,
  so `msdeconv` output (and TopDIA's default output) keeps the earlier block
  with the DP's ranking window and harmonic check lines.
- **Removed options** (unrecognised now): `-F`/`--ms2-flashdeconv`, the
  hidden `--ms1-greedy-deconv` / `--ms2-greedy-deconv`, and the short-lived
  `--deconv-ms1-msdeconv` / `--deconv-ms2-msdeconv` of commit b3d71867.
  The letter `-F` is free in topfd.
- **GUI wiring**: `topfddialog.ui` and `topdiadialog.ui` have
  `ms1DeconvComboBox` and `ms2DeconvComboBox`, whose item texts are the
  option values (as `activationComboBox`'s are), set from the dialog's
  default `TopfdPara` (TopFD: a default-constructed one; TopDIA:
  `Argument::getTopfdParaPtrForTopdia`) and read back with
  `setMs1DeconvMethod(currentText())`; `gui/util/command.cpp` always emits
  `--deconv-ms1 <m> --deconv-ms2 <m>`. In topfd_gui they are the third
  column's rows 4 and 5 (the window grew by 40 pixels, and the former "Use
  FLASHDeconv for MS/MS" checkbox is gone); in topdia_gui the middle
  column's rows 3 and 4, above the MS1 and MS2 ECScore cutoffs ("Use
  single scan noise level" moved into the Additional settings box). The
  forms place widgets at fixed coordinates, so moving one means editing the
  geometries of the widgets below it. topdia_gui's window is 1103 pixels
  wide (was 965; left column rows 305, right column rows 385, deconvolution
  drop-downs at most 120; the correlation cutoff labels say "MS1 envelope"
  and "MS2 envelope", without "isotopic") so that every label and
  drop-down fits the desktop font Ubuntu Sans 11, not only Qt's offscreen
  Sans Serif 9. To check a form, build a small program on the generated
  `build/<tool>_gui_autogen/include/gui/<tool>/ui_<tool>dialog.h`, run it
  with `QT_QPA_PLATFORM=offscreen` and the font set by
  `QApplication::setFont`, and compare each `QLabel`'s
  `sizeHint().width()` with its `width()` (and `grab().save()` a
  screenshot).

## FLASHDeconv MS/MS deconvolution (`topfd/deconv/flashdeconv_ms2_process`)

`FlashDeconvMs2Process` is a second `DeconvMs2Process`: the same mzML group
reading, one task per MS/MS spectrum on the thread pool, per-thread msalign
files merged at the end and buffered SQLite records, but the spectrum is
deconvoluted with the FLASHDeconv algorithm (OpenMS `SpectralDeconvolution`,
`PeakGroup`, `PeakGroupScoring`, `FLASHHelperClasses`) instead of
`DeconvSingleSp`. `--deconv-ms2 flashdeconv` selects it in topfd and topdia
(`TopfdPara::isUseFlashDeconv`, `topfd_process.cpp`, `topdia_process.cpp`;
previous section). Facts to keep straight:

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

## Greedy MS1 deconvolution (`topfd/deconv/greedy_ms1_deconv`, TopFD's default)

`deconv_ms1_process::deconvMsOneGreedy` is TopFD's MS1 deconvolution
(`TopfdPara::isMs1GreedyDeconv`, true by default; the parameter block then
prints "MS1 deconvolution: Greedy"). It replaces the precursor-window step
and `DeconvSingleSp` (candidate envelopes, filters, MS-Deconv's DP) of
`deconvMsOne` for MS1 spectra with `greedy_ms1_deconv::deconv`, a
THRASH-like subtractive search in which a peak's charge is decided with all
the charge states of its mass; the output (msalign lines per charge state,
EnvCNN scores, SQLite rows) is written by the same `writeMsOne` as
`deconvMsOne`'s. `deconvMsOne`, TopFD's default up to main 89df2a8c, runs
with `--deconv-ms1 msdeconv`, TopDIA's default (see "Deconvolution
methods" above for the option, the defaults and the GUIs). Facts to keep
straight:

- **One spectrum** (`GreedySearch::run`): the most intense peak of the
  residual spectrum (residual >= noise x MS1 S/N, the thresholds of
  `DeconvData::initMinInte`) is the seed; for every charge z with a peak one
  isotope spacing from the seed (pre-check), the EnvBase distribution's apex
  is aligned with the seed, scaled to the residual intensities of its core
  isotopes (>= 0.5 x apex, least squares), and extended to the same mass at
  z-1, z-2, ... and z+1, z+2, ... while those envelopes fit (isotope cosine
  >= 0.6, the apex matched, `min_match_peak_num_`; one failing charge in a
  row allowed). The charge whose series has the highest score wins; an
  envelope's score is the explained intensity (sum of min(observed,
  theoretical)) minus the theoretical intensity predicted where there is no
  peak, times the squared isotope cosine (`score_mode` 1). Then the
  monoisotopic offset (+-2 isotopes) is chosen by the cosine of the
  averagine distribution with the isotope intensities summed over the
  series' charges (`mono_by` 1, FLASHDeconv's
  `getIsotopeCosineAndIsoOffset`), and the mass's envelopes at all its
  charges are subtracted: a matched peak above 1.5x its theoretical
  intensity keeps the difference, others are used up, and the tail
  isotopes up to where the observed series stops or rises again are removed
  whole (`dyn_bound`).
- **Guards against artefacts**: an envelope within 3.5 Da of a mass already
  found at the same charge is not used (`dup_window`, as main's 12 Da
  `filterEnvByMz`; without it about 20% of a dense spectrum's masses have
  a same-charge near-duplicate, mostly from peaks the earlier envelope never
  matched rather than from its residual, and the feature detection merges
  them, so the identifications do not change); a seed charge z whose
  isotopes off the grid of every 2nd/3rd isotope are much weaker than those
  on it is rejected as a harmonic (`harm_ratio` 0.4); a mass found at one
  charge only needs an isotope cosine of 0.8 (`single_min_cos`); envelopes
  with an EnvCNN score below 0.2 are not written (`min_envcnn`; a third of
  the envelopes, which only seeded extra MS1 features).
- **Precursor windows without an envelope** (postprocessing step,
  `addPrecWinEnvs`, switch `prec_win_post`, on by default; toppic_test
  experiment 113): `deconvMsOneGreedy` passes the precursor windows of the
  MS/MS spectra that follow the MS1 spectrum (`PrecWin`, not when
  `prec_win` runs the precursor-window step first). After the EnvCNN
  filter, a window in which no reported envelope has its reference (apex)
  m/z, the test the MS1 feature detection uses to seed its search for
  windows without a feature (`env_coll_detect`, `SeedEnv::getReferMz`),
  gets the envelope of its most intense peak (`GreedySearch::fitWindow`,
  on the original intensities): one peak in the window is a charge 1
  peak; with more, every charge is tried as for a seed but without the
  acceptance thresholds (the harmonic test stays) and the best series
  score wins; the monoisotopic offset is kept only when the apex stays in
  the window. One envelope per window (the window's charge), with its
  EnvCNN score but not filtered by it; a window without a peak gets
  nothing. `TOPFD_GREEDY_PREC_WIN_LOG=<file>` writes one line per window.
  On 480_F3_02 (toppic_test experiment 113, report 0113) 80 of 4,323
  windows had no envelope (40 deconvoluted at charge 1); 20 more MS/MS
  spectra got a precursor (4,318 against 4,298), none identified (their
  MS/MS spectra have a median of 4 fragment masses), and one harmonic
  precursor was corrected (3,487 Da at charge 3 -> 6,974 Da at charge 6, a
  better PTMS proteoform). PrSMs 3,331 and 3,338 against main's 3,337 and
  3,338 (the run-to-run variation); same TopFD time. A new charge 1 feature
  can cover a later window at the same m/z, so that the feature search for
  windows without a feature skips it (scans 5806 and 5815 lost a 10,943 Da
  precursor this way; neither is identified in any run).
- **Research switches**: the defaults are in `GreedyPara` (cpp); the
  environment variable `TOPFD_GREEDY_PARA` ("name=value,...") overrides
  them, so that variants can be run without a rebuild (toppic_test
  experiment 104 lists the runs). `TOPFD_GREEDY_DEBUG_SCAN=<scan>` writes
  the decisions of one MS1 scan to `greedy_scan<scan>.txt` in the working
  directory.
- **Cost**: the search itself is cheap (about 70 ms per 480_F3_02 MS1
  spectrum); almost all of the extra MS1 time against main is EnvCNN on
  the envelopes before the `min_envcnn` filter (about twice main's
  number).
- **Results** (toppic_test experiment 104, report 0104): on 480_F3_02
  with main's TopPIC in the benchmark setting, 3,181-3,182 PrSMs, 1,887
  proteoforms and 425 proteins against main's 3,141-3,153 / 1,870 / 420;
  the gain is corrected precursor isotopes and charges (and fewer decoys
  on the spectra whose precursor changed). TopFD takes about the same time
  (334-336 s against 285-333 s) and memory. With both greedy searches
  (toppic_test experiment 109, rounds 2 and 3): 3,316-3,318 / 1,981 / 433.
  The search gives fewer MS/MS spectra a precursor on some data: 612 of
  1,381 against MS-Deconv's 855 on st_1.mzML (MS/MS-matched features 894
  against 1,158), 4,289 against 4,303 on 480_F3_02.
- **`--deconv-ms1 msdeconv` writes the earlier default's output**: on st_1
  all 1,418 MS1 scans' masses and charges, both feature files and the
  parameter block equal main 89df2a8c's (the last main with the DP as the
  default), the MS/MS spectra up to the EnvCNN jitter of two runs of one
  binary. The greedy search itself is not exactly repeatable (about 5% of
  st_1's MS1 scans differ between two runs).
- **Not done**: the precursor-window step of `deconvMsOne` (switch
  `prec_win` runs it first; it did not change the identifications) and
  FAIMS/`-T` specifics (the MS1 path is the same). The MS/MS spectra have
  their own greedy search (next section).

## Greedy MS/MS deconvolution (`topfd/deconv/greedy_ms2_deconv`, TopFD's default)

`greedy_ms2_deconv::deconv` is TopFD's MS/MS deconvolution
(`TopfdPara::isMs2GreedyDeconv`, true by default; `deconv_ms2_process::
deconvMsTwo` and `topfd -T` call it; the parameter block then prints
"MS/MS deconvolution: Greedy" instead of the DP's ranking window and
harmonic check lines). `DeconvSingleSp`, TopFD's default up to main
89df2a8c, runs with `--deconv-ms2 msdeconv`, TopDIA's default (see
"Deconvolution methods" above). `-v`, `-g` and `-n` apply to both methods
(`greedy_ms2_deconv::deconv` ends with the same filters). It started as a
copy of the greedy MS1 search (previous section) and keeps its structure;
the output (msalign lines, EnvCNN scores, the `-F`-relevant score, SQLite
rows) is written as `DeconvSingleSp`'s. Facts to keep straight:

- **Bounds**: fragment masses below the precursors' largest mass and
  charges up to their largest charge (`deconvMsTwo`'s `max_mass` /
  `max_charge`), the MS/MS signal-to-noise ratio for the seed threshold.
- **What differs from the MS1 search** (each a `GreedyMs2Para` switch,
  toppic_test experiment 106):
    - every charge competes with the best of its isotope offsets
      (`charge_offsets`): below about 5,000 Da the averagine apex and the
      monoisotopic peak are of similar height, and the apex alignment of the
      right charge was often one isotope off and lost to a wrong charge;
    - the offset is chosen on the seed charge's envelope alone
      (`offset_by_seed`), and a neighbouring charge joins the series only
      when its own best offset agrees (`nb_consistent`); a neighbouring
      charge was often the envelope of another mass one isotope away;
    - the search matches peaks within 5 ppm (`ppm_tol`), not TopFD's 0.02
      m/z: on 480_F3_02 the peaks of confident envelopes lie within 0.6 ppm
      (median), and the wide window let high-charge envelopes collect the
      peaks of other species; a peak only within 0.02 m/z (two envelopes'
      peaks merged into one centroid) is neither matched nor counted as
      missing (`near_no_pen`). The reported envelopes are matched with
      TopFD's tolerance, as EnvCNN expects;
    - the same-charge duplicate test only blocks masses an integer number
      of isotope spacings from an earlier one (`dup_iso`): distinct
      masses of one charge can lie closer than 2 Da;
    - no harmonic test and a single-charge cosine of 0.7 (most fragments
      appear at one charge);
    - the charge pre-check also accepts a neighbour peak at the averagine
      distribution's own isotope spacing (`precheck` 2): the M+1 peak of
      averagine is 1.00286 Da from the monoisotopic peak at every mass,
      while the isotope mass constant (`mass_constant::getIsotopeMass`,
      1.00228 Da, the apex spacing of human averagine envelopes above about
      5 kDa; 1.00235 before toppic_test experiment 111) fits the apex
      (toppic_test experiments 107 and 108); at m/z 200 and charge 1 the
      difference is 2.9 ppm, more than half of the 5 ppm window.
- **Research switches**: `TOPFD_GREEDY_MS2_PARA` ("name=value,...")
  overrides the defaults in `GreedyMs2Para`; `TOPFD_GREEDY_MS2_DEBUG_SCAN=
  <scan>` writes one scan's decisions to `greedy_ms2_scan<scan>.txt`.
  Switches kept off by default are listed there with their experiment 106
  results (`refit_rounds`, `alt_ratio`, `mono_by`, `score_mode` 2, ...).
- **Results** (toppic_test experiment 106, report 0106): on 480_F3_02
  with main's TopPIC in the benchmark setting, 3,296-3,298 PrSMs, 1,961-
  1,964 proteoforms and 427 proteins against main's 3,141 / 1,870 / 420 on
  the same MS1 output (+199 to +228 targets at E <= 1e-5 to 0.01, no more
  decoys); with the greedy MS1 search too 3,315 / 1,980 / 433. The gain is
  15% more b/y fragments of the identified proteoforms at main's precision
  (mostly 1,000-5,000 Da: right charges and isotopes, and envelopes the DP
  did not report). The MS/MS stage takes about 100 s against main's 135 s
  at 16 threads.
- **`--deconv-ms2 msdeconv` writes the earlier default's output**: on
  st_1, with `--deconv-ms1 msdeconv` too, the MS/MS spectra equal main
  89df2a8c's default output up to the EnvCNN jitter, and `topfd -T` on the
  myoglobin peak list writes main's msalign byte for byte (681 masses; the
  greedy search 2,476). With the MS1 deconvolution fixed, the greedy MS/MS
  search is repeatable up to that jitter. `--deconv-ms2 flashdeconv` writes
  main 89df2a8c's `-F` MS/MS spectra unchanged (st_1, all 1,381).

## Precursor isotope correction from fragment pairs (`topfd/deconv/prec_pair_refine`, `--skip-ms2-prec-pair-refine`)

On by default in topfd; the hidden `--skip-ms2-prec-pair-refine` turns it
off (branch `high_res_prec_pair_refine`, toppic_test experiment 116, report
0116; tested as the hidden `--ms2-prec-pair-refine`, which the merge into
main replaced; `--skip-ms2-prec-pair-refine` was visible in `-h` until
2026-10-08, when it was moved to the hidden options).
b + y = the proteoform mass, so fragment pair sums show the mass of the
molecule that was fragmented. Facts to keep straight:

- **Two functions.** `prec_pair_refine::refine` is the low-resolution
  method (`-R`, branch `low_res_greedy_ms1` f902aa09, copied unchanged; main
  has no `-R`, so nothing calls it here; a later merge of a low-resolution
  branch has to rename its `--ms2-prec-pair-refine` too). `correctIsotopes`
  is the high-resolution one that TopFD runs, in `deconvMsTwo` after the
  precursors are built from the features.
- **Where it runs.** `TopfdPara::ms2_prec_pair_refine_` is true by default.
  It runs in `DeconvMs2Process` only (greedy and msdeconv MS/MS
  deconvolution, not FLASHDeconv or `-T`) and only for MS/MS spectra with an
  MS1 spectrum: without one (`-o`) the MS1 check cannot run, and with the
  pairs alone the rule moved 47 right precursors of the known-mass spectra
  away and fixed 18 (report 0116 section 3).
  TopDIA's argument parser turns it off (DIA MS/MS scans hold the fragments
  of every co-isolated precursor; not tested). The parameter line
  `Precursor isotope correction from fragment pairs: Yes` is printed only
  when the correction runs, so the other output files keep the parameter
  block of earlier versions.
- **The rule (`correctIsotopes`).** For every precursor of a spectrum (TopPIC
  searches only the first), it counts the pair sums of the fragment masses
  with EnvCNN >= 0.4 within 10 ppm of P - 1, P and P + 1 isotope (1.00228 Da).
  The precursor moves to the best offset when it has >= 3 pairs, more than at P
  and >= 2 x the count at P. It does not move when the averagine fit of the
  precursor in the group's MS1 spectrum prefers P by 0.05 or more (cosines on
  one grid of isotope positions, +-2 isotopes around the envelope, summed
  over the charges z-1..z+1).
    - The moved precursor keeps its charge, feature id and intensity; its m/z
      moves by j x 1.00228 / z.
    - `DeconvMs2Process::process` passes the group's MS1 spectrum to the
      MS/MS task only when the correction is on.
- **Research switches** in `TOPFD_PREC_REFINE_PARA`: `iso_range` (1),
  `iso_tol_ppm` (10), `iso_min_envcnn` (0.4), `iso_min_pairs` (3),
  `iso_ratio` (2), `iso_ms1_check` (1), `iso_ms1_margin` (0.05),
  `iso_ms1_tol_ppm` (10), `iso_ms1_min_rel` (0.05). The low-resolution
  `window`, `tol_ppm`, `min_pairs` and `promote` are unchanged.
- **Results on 480_F3_02** (`-s 2`, oxidation only):
    - It moves 96 first precursors. At 1% FDR TopPIC finds 3,341 PrSMs /
      2,003 proteoforms / 438 proteins, against main's 3,337 / 1,997 / 438.
    - The moved spectra get better PrSMs: of 86 identified by either run,
      56 match more fragments and 10 fewer.
    - 12 of them are Asn-Gly deamidations (HMGN1, parathymosin, hnRNP A1-like
      2) that main reports without the deamidation. A deamidated proteoform is
      0.018 Da short of one isotope heavier, and MS1 cannot separate its
      envelope from the unmodified one.
- **Why the rule is conservative.** TopFD's fragment masses carry isotope
  errors, mostly one isotope light: they are right in 98% of the matches
  below 1 kDa but 78% at 3-4 kDa and 55% at 6-10 kDa. A pair sum has the
  isotope error of its larger fragment. With pairs alone (no MS1 check,
  moves of +-2) 240 precursors move and many PrSMs get worse.
- **TopFD's own precursor errors** are rare: 97 of 1,866 known-mass spectra,
  half of them histone H2B at 13.8 kDa (one isotope light). The pairs cannot
  reach those: complementary pairs are missing or flat there.

## PTCR spectra (`topfd/deconv/ptcr_ms2_deconv`)

Branch `ptcr_ms2_deconv` (toppic_test experiment 136, report 0136; merged
into main on 2026-10-07 as 5f19b7cb: with the branch's float EnvCNN / ECScore
models the merged topfd writes the branch's PTCR, fragment PTCR and LC-MS
outputs, and main's quantized models change only EnvCNN scores there, no
mass line; 480_F3_02 within the variation of two runs of main faa8e9c3;
report 0136 section 9). A PTCR spectrum is an MS/MS spectrum of proton transfer charge reduction of the
isolated ions without dissociation: no fragments, only the isolated intact
species at lower charges. Facts to keep straight:

- **Reading.** `PwMsReader::parseActivation` decides on the activation
  terms of all precursors of a spectrum (an MSn scan has one per stage;
  msconvert writes the first stage first). `getActivationName` names the
  dissociation (proton transfer is none); `hasPtcrTerm` finds MS:1003249
  (proton transfer charge reduction) or MS:1003248 (proton transfer
  reaction). A proton transfer term without a dissociation: a PTCR
  spectrum (`MsHeader::isPtcr`, an activation object of its own,
  `getPtcrActivationPtr`, placeholder ion types B / Y). `PTCR` is
  deliberately not in `activation_base.xml`, the `-a` lists or the GUIs: no
  search reads these spectra. A dissociation (of the file or of `-a`) and a
  proton transfer term: a fragment PTCR spectrum (`isFragPtcr`, below), of
  the dissociation's activation, when the proton transfer stage comes after
  the last dissociation stage or in its stage; so `-a HCD` on a PTCR scan
  now makes a fragment PTCR spectrum, not main's fragment spectrum. **The
  order of the stages** decides (toppic_test experiment 139):
  `getScanStages` lists them in acquisition order (the filter string's
  stages, each with the file's precursor of the same target, 0.01 m/z;
  without a filter string the precursors in file order, msconvert's stage
  order), and when a dissociation stage comes after the last proton
  transfer stage (`getReducedPrecStage`) the spectrum is the dissociation
  of a charge-reduced precursor (`isReducedPrec`, below), not a fragment
  PTCR spectrum. Main 6eb40e97 read the MS3 scans "X@ptr5 Y@hcd80" of
  PRIDE PXD074730 as fragment PTCR spectra with the window of the proton
  transfer stage and, without a charge, `-o`'s placeholder precursor. The hidden `topfd
  --ptcr` sets the activation to `PTCR` for every MS/MS spectrum
  (`TopfdPara::isPtcrForced`, which the activation check accepts), for
  files whose converter wrote no proton transfer term. `parsePrecursor`
  keeps the selected ion's m/z and charge of the first precursor
  (`MsHeader::getSelectedMz` / `getSelectedCharge`, -1 / 0 without) and the
  isolation window of the proton transfer stage (`getPtcrWinBegin` /
  `getPtcrWinEnd`: the last precursor with a proton transfer term, else the
  precursor window).
- **MSn scans** (n > 2) with a proton transfer stage are MS/MS spectra
  (`MzmlMsGroupReader`'s `isMsTwo`, in `-o` and with MS1 spectra, and the
  count of `readProfile`): grouped with the last MS1 scan before them
  (`parsePrecursor` ignores their precursor reference, the MS(n-1) scan).
  Other MSn scans (a dissociation of fragments, whose precursor is not the
  window's proteoform) are skipped, as main skipped every MSn scan.
  **Missing stages** (`getMissingStages`): ThermoRawFileParser 2.0 writes
  one precursor for the MS3 scan "ms3 1402.5000@cid60.00
  1100.0000@ptr10.00", the proton transfer stage (window 1100 +- 900, the
  PTR term, the scan's charge state), and loses the CID. When a scan has
  fewer precursors than stages and its filter string (MS:1000512) lists
  them, the stages' terms (cid, hcd, etd, ptr) join the activation terms,
  and when the file's precursor is not the first stage, the precursor
  window and the selected ion m/z are the first stage's target with the
  default window (`-w`); the charge stays the file's. The same holds for
  **MS2 scans without a precursorList** (toppic_test experiment 139):
  ThermoRawFileParser writes none for a scan whose precursor scan is not in
  the file (a scan cut out of a longer run, as all MS2 and MS3 files of
  PXD074730), and main stopped on them ("Unable to read the activation
  method"); the filter string's stage now gives the target with the
  default window, the selected m/z (no charge) and the activation.
- **Where it runs.** `deconv_ms2_process::deconvMsTwo` sends a PTCR
  spectrum to `ptcr_ms2_deconv::deconv` and buffers the result
  (`PtcrBuffer`, under a lock); `DeconvMs2Process::process` writes
  `<base>_ptcr.msalign` (msalign format, `ACTIVATION=PTCR`, one line per
  species: mass, intensity, highest charge, best EnvCNN score; `-A` average
  masses) and `<base>_ptcr.tsv` (one line per species with its charges, the
  intensity and EnvCNN score per charge, the isolation charges, the
  precursor flag) after the pool, sorted by spectrum id, only when there are
  PTCR spectra. They are not in `_ms2.msalign` or the SQLite database. In a
  run with MS1 spectra, a species whose mass matches an MS1 feature of the
  spectrum's window (`feature_ppm` 10 ppm, also one isotope apart) gets its
  id and mass (the TSV columns `Feature ID`, `Feature mass`), and the
  precursor species' feature gives the precursor its id and intensity.
  Not in `FlashDeconvMs2Process` (`--deconv-ms2 flashdeconv`) or `-T`.
- **The search** is the greedy MS1 search through
  `greedy_ms1_deconv::deconvSeries`, a second entry point of
  `greedy_ms1_deconv` (no precursor-window step) that takes a seed S/N
  (the MS/MS ratio `-s`) and a `ChargeBound` (`fitEnv` rejects an excluded
  (mass, charge): seeds, series and offsets) and returns the envelopes
  grouped by mass (EnvCNN scored, below 0.2 dropped, intensities shared
  among the kept ones). The PTCR bound: charge z of mass M only when some
  z0 >= z puts the envelope of M (its isotopes >= 10% of the apex in
  TopFD's averagine table) in the isolation window, since proton transfer
  only lowers charges. Without an isolation window there is no bound.
- **A species** is one mass's envelopes merged: the summed intensity, the
  charges, the mass refined from the matched isotope peaks of all charges
  (the intensity-weighted median of (observed - theoretical m/z) x z; the
  envelope's own mass rests on the seed peak). **The precursor** is a
  species with an envelope at the selected ion's charge whose monoisotopic
  m/z at that charge is the selected ion's m/z (`sel_ppm`, 10 ppm), else
  the most intense species with an envelope at that charge whose m/z
  range, from its monoisotopic peak to its last isotope of >= 10% of the
  apex, covers the selected ion's m/z, else the most intense one that could
  have been isolated at that charge, else the most intense one, else the
  selected ion itself. Thermo converters write the instrument's
  monoisotopic m/z (the "Monoisotopic M/Z" trailer) as the selected ion
  when it has one, else the isolation target. The cover test began at the
  first isotope of >= 10% (toppic_test experiment 136), at 17 kDa several
  isotopes above the monoisotopic peak, so a monoisotopic selected m/z lay
  in the envelope of the 16 Da lighter species only: on PXD074730's
  Lysozyme_MS2_PTCR (selected 941.9354, 18+: 16,936.70 Da) main 6eb40e97
  took the 16,920.73 Da species (6.9 M) instead of the 16,936.69 Da one
  (80.7 M; experiment 139).
- **Fragment PTCR spectra** (`isFragPtcr`: a dissociation, then proton
  transfer charge reduction of the fragments) are deconvoluted as other
  MS/MS spectra of their activation (`deconvMsTwo`, `--deconv-ms2`) and
  written to `_ms2.msalign` and the SQLite database, so TopPIC searches
  them with the dissociation's ion types. The selected ion bounds the
  fragments (`ptcr_ms2_deconv::getSelectedIonPrec`: its charge, and its
  mass at the upper end of the precursor window; not capped by `-c` /
  `-m`), with MS1 features the larger of these and the features' bounds (a
  feature of the window need not be the isolated ion: the AAV2 MS1 scans
  are not isotope-resolved and their features are noise, which as the only
  bound left 470 masses per spectrum instead of 4,100-5,500). Without MS1
  features of the window the selected ion (its m/z and charge as the file
  gives them, not monoisotopic) is the precursor. The console counts them
  ("Fragment PTCR spectra: N"). Switch `frag`: 1 this (default), 0 as other
  MS/MS spectra (`-o`: the dummy precursor of `-c` / `-m`), 2
  `deconvFragments`, the charge-series search of the intact spectra with
  the PTCR bound, z0 <= the selected charge and masses below the
  precursor's (also its surviving precursor ladder as the precursor), one
  mass per fragment species, `<base>_ptcr_fragment.tsv`. The MS/MS search
  matched more theoretical fragments of AAV2's capsid proteins than the
  series search (below): with the wide PTCR windows of these data
  (1100 +- 900, 650 +- 255 m/z) the band bound excludes almost nothing
  above 1 kDa, and the MS/MS search's fragment-specific choices (5 ppm,
  per-charge isotope offsets, no harmonic test, the single-charge cosine
  0.7) find more low-charge fragments.
- **The dissociation of a charge-reduced precursor** (`isReducedPrec`;
  toppic_test experiment 139): an MSn scan whose dissociation stage comes
  after its proton transfer stage, e.g. the MS3 scan "X@ptr5 Y@hcd80" of
  the cycle MS1, PTCR MS2 of X, HCD MS3 of one charge-reduced ion at Y
  (PXD074730, iodoTMT quantification of intact proteins: the MS3 scans
  hold the reporter ions). It is read as an MS/MS spectrum of the
  dissociation's activation (`isMsTwo`, `isReadMsTwo`) and written to
  `_ms2.msalign` (`LEVEL=3`) and the SQLite database; its precursor window
  is the first stage's, the MS1 window of the isolated ions, matched to
  the MS1 features as for every MS/MS spectrum. The header keeps the
  dissociation stage's target and window (`getReducedPrecTargetMz`,
  `getReducedPrecWinBegin` / `End`: the file's precursor of that target,
  else the filter string's target with the default window `-w`) and the
  scan its last precursor references (`getParentScan`). Its precursor
  (`deconv_ms2_process::deconvMsTwo`, steps 0b and 1c): (a) from the
  parent PTCR spectrum (`findPtcrParent`: the referenced scan, else the
  last PTCR spectrum read before it with the same isolation target, 0.01
  m/z, among the last 64), deconvoluted again in the MS3's task
  (`ptcr_ms2_deconv::getReducedPrec`; `getPtcrSpec` reads the parent's
  header and does not change it, so both tasks can run at once): the most
  intense species with an envelope at one of its charges, else at any
  charge up to its isolation charge, whose m/z range (monoisotopic peak to
  the last isotope >= 10%) contains the dissociation stage's target, at
  that charge, with its intensity and MS1 feature; (b) without a parent,
  the MS1 features of the window with an envelope at a charge up to theirs
  containing that target (`getReducedPrecFeatures`), all of them when none
  has one; (c) without MS1 features the selected ion with its charge
  (`getSelectedIonPrec`); (d) else `-o`'s placeholder. The precursor's mass
  and charge bound the fragments. The console reports "Dissociation
  spectra of charge-reduced precursors: N, M with the precursor from a
  PTCR spectrum". The MS3 scan does not tell the charge of the ion at Y;
  the parent's charge ladders do. The precursor mass is the parent
  species' own, so it shares the species' isotope error (below).
- **Research switches** in `TOPFD_PTCR_PARA` ("name=value,..."): `band`
  (1), `win_margin` (0 m/z, added to both window ends), `iso_min_rel`
  (0.1), `min_envcnn` (0.2), `sn_ratio` (0 = `-s`), `refine` (1), `frag`
  (1, above), `frag_min_envcnn` (0, frag 2), `frag_prec_margin` (20 Da,
  frag 2: species this close to the precursor mass or above are not
  fragments), `feature_ppm` (10), `sel_ppm` (10, the monoisotopic
  selected ion). The greedy search's own switches
  (`TOPFD_GREEDY_PARA`) apply too.
  `TOPFD_PTCR_DEBUG_SCAN=<scan>` writes one scan's decisions to
  `ptcr_scan<scan>.txt`.
- **`-o` lost the first MS/MS spectrum of every file without MS1
  spectra**: `MzmlMsGroupReader`'s constructor reads one spectrum to check
  that the file is not empty and rewound only when MS1 spectra are present
  (fixed in this branch). TopDIA's MS/MS step runs in that mode on files
  that start with an MS1 spectrum, which the skip consumed harmlessly.
- **Results** (toppic_test experiment 136; data/ptcr_ms2, two single PTCR
  scans of an Orbitrap Eclipse, isolated at 16+ and 17+ with +-5 m/z, of a
  16.94 kDa proteoform and its co-isolated variants; ThermoRawFileParser
  drops the precursor of these scans, cut out of longer runs, and the
  experiment's add_precursor.py writes it back): the dominant proteoform at
  16,936.657 / 16,936.663 Da in the two scans (0.3 ppm; 0.9 ppm without the
  refinement), 11+-16+ and 11+-17+, as their precursor; 77 / 95 species,
  no 1-isotope twins; species with 75% / 72% of the intensity found in both
  scans at the same mass (median 0.47 ppm apart), 17% / 16% at an isotope
  shift. The co-isolated variants lie a near-whole number of Daltons apart,
  so their isotope peaks merge into one centroid (5 ppm at every charge),
  and a minor species' monoisotopic mass is not determined to the isotope
  (a joint least-squares fit of all species does not decide it either).
  Without the bound 81 / 102 more species (9% / 10% of the intensity),
  2 of them at the same mass in the other scan. One species whose envelope
  ends 0.8 and 0.5 m/z below the two windows is reproducible (16,835.55 Da,
  1% of the intensity): the window edges are not sharp, but margins of
  0.5-3 m/z added mostly irreproducible species, so the margin stays 0. main 38ba67dc writes 420 /
  487 envelope lines (one per charge state) with a dummy precursor (`-o`:
  charge 30, mass 50,000 Da), 21 / 33 masses above 10 kDa with a
  1-isotope twin, 72% / 77% of the lines (15% / 19% of the intensity) of
  masses outside every band, charges up to 30. On 480_F3_02 (no PTCR
  spectra) the branch's TopFD output equals main's within the variation of
  two main runs (MS/MS spectra identical; 702 against 727 MS1 spectra with
  a differing mass; same time and memory).
- **Results of fragment PTCR and LC-MS runs** (toppic_test experiment 136,
  amendment 1; MassIVE MSV000095045, AAV2 capsid on an Orbitrap Ascend: 36
  MS3 scans "ms3 1288|1300|1402.5@cid60 1100@ptr10" and 36 "825|830@etd3
  650@ptr10|8", matched against the b/y or c/z-dot ions of VP1, VP2 and VP3
  at 10 ppm with reversed-sequence decoys): the MS/MS search (frag 1) 1,125
  CID and 334 ETD target matches above the decoys' (sums over the
  spectra), the series search (frag 2) 937 and 206, frag 0 with `-o`'s
  bounds 1,113 / 320 and with `-c 64 -m 90000` 1,122 / 335. Main reads no
  MS3 scan (0 MS/MS scans). ThermoRawFileParser's own files (first stage
  from the filter string) give the masses of the files with a precursor per
  stage in 69 of 72 spectra (3 ETD spectra differ near the bound: the
  default window of 3 m/z against +-2.5). A synthetic LC-MS run (the
  AAV2_MS1 run's 2,882 MS1 scans with the 72 MS3 scans and 3 MS1 scans
  copied as intact PTCR scans): every MS3 scan in the group of the MS1 scan
  before it, the masses of the `-o` runs in all 72; 62 precursors from
  features (noise: 43 of charge 1), 10 the selected ion; an intact copy's
  species linked to its window's feature. The intact PTCR outputs of
  data/ptcr_ms2 are unchanged (but the two TSV columns); 480_F3_02 within
  the variation of two main runs.
- **Results on PRIDE PXD074730** (toppic_test experiment 139, report 0139;
  Kaulich and Tholey 2026, iodoTMT-labeled intact lysozyme, BSA and
  beta-lactoglobulin by direct infusion on an Orbitrap Eclipse: MS1, PTCR
  MS2 of X, HCD MS3 of one charge-reduced ion at Y with the reporter ions;
  19 single scans cut out of runs, ThermoRawFileParser 2.0.0-dev, 18 without
  a precursorList). Main 6eb40e97 stopped on the 11 MS2 files without a
  precursorList, took the 16,920.73 Da species as the precursor of
  `Lysozyme_MS2_PTCR` (selected 941.9354 at 18+ = 16,936.69 Da, the
  strongest species), and read the 7 MS3 scans as fragment PTCR spectra
  with the placeholder precursor. The branch reads all 19; with each PTCR
  scan and its MS3 scans in one file, the 6 MS3 precursors come from the
  PTCR species at Y (lysozyme 16,936.66-16,936.69 Da at 15+ to 18+, within
  0.6 ppm of the instrument's monoisotopic m/z of Y in the MS3 scan's
  trailer; beta-lactoglobulin 19,915.54 Da at 18+, one isotope above the
  instrument's mass and its MS1 feature), also from the stage precursors
  alone (no filter string) and with a spectrumRef; with MS1 scans (one
  run's 396 MS1 scans with the PTCR scan and two MS3 scans inserted) each
  MS3 gets its own proteoform and feature, where main gave both scans both
  features. Fragment PTCR on AAV2 and the 4 MS1-only files unchanged
  (masses; EnvCNN scores by up to 0.02); on 480_F3_02 the reading fixes
  alone (7e5b5a64) are within the variation of two main runs, and with the
  EnvCNN input (above) toppic's identifications equal main's.

## MS3 spectral deconvolution (`topfd/deconv/deconv_ms3_process`, `--disable-ms3`)

Branch `ms3_deconv` (toppic_test experiment 146, report 0146; from main
205dc03c, merged into main as ceab7362 on 2026-10-08). An MS3 scan is the
dissociation of an MS2 fragment ion (a
tryptic digest on an Orbitrap Fusion Lumos: CID MS2 of the peptide, CID
MS3 of up to three of its fragments). TopFD read only the MS3 scans with a
proton transfer stage (PTCR, above); the others were skipped. Since this
branch they are deconvoluted by a stage after the MS/MS deconvolution,
`DeconvMs3Process`, on by default when the file has such scans
(`MzmlProfile::getMs3Cnt`, counted by `PwMsReader::readProfile`; the
parameter block prints "Number of MS3 scans" and "MS3 spectral
deconvolution: on / off"; hidden `topfd --disable-ms3`; not run for
FAIMS data). Facts to keep straight:

- **Reading.** The stage reads the mzML file again with `PwMsReader`
  (which reads every scan; the MS/MS stage's `MzmlMsGroupReader` is what
  skips MS3 scans) and keeps the level-3 scans that are not PTCR,
  fragment PTCR or charge-reduced-precursor spectra. The parent MS2 scan
  of an MSn scan (`MsHeader::getParentScan`) is the `spectrumRef` of the
  precursor element whose isolation target is the **last** stage of the
  filter string (`getFilterStages`): ThermoRawFileParser 2.0 writes the
  MS3 stage first and the MS1 stage second, msconvert the MS1 stage
  first; without a filter string the last precursor (msconvert's order).
  The MS1 scan is the cycle's (`MS_ONE_SCAN`, `prev_ms1_scan_id_`). The
  isolation target and window are the first precursor element's
  (`parsePrecursor`), which is the MS3 stage for ThermoRawFileParser
  files - with msconvert files it would be the MS1 stage (not tested).
  The activation is the stages' (CID).
- **The precursor** of an MS3 spectrum is an MS2 fragment: the species
  of the parent MS2 scan's deconvolution (`<base>_ms2.msalign`, read
  back by scan) whose isotopic envelope (the monoisotopic m/z and the
  `1 + mass / 800` isotope peaks above it) has a peak in the MS3
  isolation window, the most intense when several (`findPrecSpecies`);
  its mass, charge, m/z and intensity are the precursor
  (`PRECURSOR_FEATURE_ID` -1: a fragment has no MS1 feature). Without
  one - in the test data the MS3 scans of the noise MS2 scans at the
  gradient's ends, whose MS2 precursor has the dummy charge 1 - the
  selected ion m/z of the file at every charge from 1 to the MS2
  precursor's charge, written as multiple precursors (the file's charge
  state of an MS3 precursor is not usable: 10 for a third of the scans).
  On the test data the chosen species is a b or y ion of the identified
  parent peptide in two thirds of the cases and one isotope off (the
  MS2 deconvolution's monoisotopic error: -1 Da for 2-6 kDa fragments at
  charge 3-6, +1 Da below 3 kDa) in a quarter, so the species' mass is
  **corrected by one isotope from the parent MS2 scan's centroid peaks**
  (`correctIsotope`; the stage keeps the MS2 scans' peaks as it reads
  them): when the peak one isotope above the reported monoisotopic m/z
  is more than `TOPFD_MS3_ISO_K` (4) times the theoretical ratio
  (`EnvBase::getEnvByMonoMass`) times the reported peak, the mass goes up
  one isotope; else when a peak one isotope below holds at least
  `TOPFD_MS3_ISO_C` (0.1) of the reported peak, it goes down one;
  `TOPFD_MS3_ISO_CORRECT=0` turns it off and the `_ms3.tsv` column
  "Isotope correction" records the shift. A cosine of the whole
  theoretical envelope against the peaks was tried first and did not
  separate the cases (experiment 146, amendment 1). TopPIC's zero-shift
  search tolerates the +-1 Da precursor error anyway.
- **The deconvolution** is the MS/MS one (`greedy_ms2_deconv::deconv`,
  or `DeconvSingleSp` at level 2 with `--deconv-ms2 msdeconv`) bounded by
  the precursor's mass (+1 Da) and charge (the MS2 precursor's in the
  fallback case), single-threaded (8,000 MS3 scans of 75 peaks take
  about 80 s).
- **Output.** `<base>_ms3.msalign` in the msalign format with `LEVEL=3`,
  `MS_ONE_ID` / `MS_ONE_SCAN`, and the new lines `MS_TWO_ID` /
  `MS_TWO_SCAN` (the parent MS2 spectrum's id in `_ms2.msalign` and its
  scan) and `MS_TWO_PRECURSOR_MZ` / `_CHARGE` / `_MASS` (the parent's
  first precursor; `MsHeader::setMsTwoId`, `setMsTwoPrec`; written by
  `MsAlignWriter` for level 3, ignored by `MsAlignReader`), the
  precursor lines and the masses; `<base>_ms3.tsv`, one row per MS3 scan
  (scan, MS2 scan and id, MS1 scan and id, retention time, the MS2
  precursor, the isolation target and window, the precursor source
  `ms2_species` / `selected_ion`, the precursor(s), peaks, masses). The
  MS1 / MS2 msalign and feature files and the SQLite database are
  unchanged (checked file for file against main). TopPIC searches an MS3
  msalign file with `-x` (no feature file) and `--disable-post-match`
  (the SQLite database has no MS3 peaks); a y-ion precursor is a
  both-truncated proteoform with no shift, a b-ion one carries -18.011
  Da as its unexpected shift.

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
  `--topmg pipeline` or `topmg-only`, max(`-s`, `--topmg-num-shift`)
  unexpected shifts, since a graph-alignment proteoform is corrected like
  any other). A corrected PrSM
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
- **Result file names carry `_post_ms2`** whenever post matching ran
  (`postBaseName`: the database's base name, so also for the `_calib` and
  `_mono` copies of `--mass-calib` and `-g`, with the `_resid` suffix of
  `--second-precursor`'s residual spectra kept: `<base>_resid_post_ms2`); the
  cleanup in `toppic_process` (`cleanToppicDir`) takes that name to find the
  intermediate files. With `-D` the names are the plain `<base>_ms2` ones.
- **Combined runs (`-c`) merge the spectra and record files, re-cluster,
  compute one FDR and report.** `TopPICProgress_multi_file` merges each
  fraction's post-matched spectra (`<frac>_post_ms2.msalign`, or
  `<frac>_ms2.msalign` with `-D`) into `<name>_ms2.msalign` with
  `msalign_frac_merge::mergeFractions`, and each fraction's
  `toppic_record_prsm` (`<frac>_post_ms2.toppic_record_prsm`, or
  `<frac>_ms2.…` with `-D`) into `<name>_ms2.toppic_raw_prsm`; `<frac>` is
  the mzML base without the `_calib` suffix of `--mass-calib` and the
  `_mono` suffix of a `-g` run (the same rule as
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
- **With `-c` the files are reported after the combined result, with its
  proteoform ids.** `TopPICProgress` skips the file's `TopPIC_report`
  when `combinedOutputName` is set, and `TopPICProgress_multi_file` defers
  the file's `cleanToppicDir` (the report needs the post-matched spectra
  and their index, and the cleanup copies the report's cutoff files to the
  XML outputs). After the combined report,
  `prsm_util::updateProteoClusterIds` (the counterpart of
  `mergePrsmFiles`) gives every PrSM of each file's `toppic_record_prsm`
  the `proteo_cluster_id` of its combined PrSM (spectrum id
  `SpPara::getMaxSpecNumPerFile()` x file index + its own, same precursor
  id, protein, start and end; a PrSM without one, not expected, gets a new
  id above the combined ones with a warning) and rewrites the file; then
  each file gets `geneSpIndex` on its spectra (`<frac>_post_ms2.msalign`,
  or the input file with `-D`), `TopPIC_report` (start time = the combined
  run's, as before) and `cleanToppicDir`. With `-C` the same runs on the
  record files of the earlier run, so their reports are rewritten. Only
  the proteoform ids change: the files' FDRs (the proteoform-level FDR of
  a PrSM is still that of its cluster in the file), `prot_cluster_id`s and
  everything else stay the file's own, so the PrSM and protein selections
  are those of a run without `-c`; the proteoform table picks the best
  PrSM per combined id among the PrSMs that pass, so where the combined
  clustering (`prsm_simple_cluster`) splits or joins the file's clusters
  (`prsm_feature_cluster`) it lists more or fewer proteoforms. TopMG's
  combined runs (`topmg_process`) are unchanged. Results (toppic_test
  experiment 122, 480_F3_02 as two TopFD replicates and as two retention
  time parts at 0-55 / 45-100 min, setting d with `-c`): every row of the
  files' record files, tables and SQLite tables has its combined PrSM's id;
  without the ids the record files, PrSM and protein tables equal main's;
  the combined outputs, `-C` reruns and a run without `-c` equal main's or
  the full run's; the proteoform tables gain 5-6 rows per file where the
  file's feature clustering joined same-protein, same-mass PrSMs of
  barely overlapping ranges (HMGN2 28-56 / 60-90) and lose 1 where it
  compares adjusted precursor masses 1.25 Da apart and the combined one
  the observed ones 1.02 Da apart (PTMA). Combined proteoforms reported
  in both files under one id: replicates 2,002 of 2,002 (main 1,658 of
  1,996), parts 404 of 404 (main 2 of 403). Time unchanged (id update
  ~1 s).
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
  reports 0023/0024). TopMG's identification is hosted behind the hidden
  `--topmg <disabled|pipeline|topmg-only>` (argument key `topmgMode`,
  default `disabled`; any other value is refused). With `pipeline`,
  `TopPIC_topmg_identify` runs after the multiple-shift stage on the spectra
  the earlier stages left unidentified and hands `toppic_topmg_raw_prsm` to
  `TopPIC_postprocess`; with `topmg-only` (which needs `-b`)
  `TopPICProgress` skips `TopPIC_identify` and
  `TopPIC_multi_ptm_search` altogether, runs `TopPIC_init` (start time,
  argument echo, base data, feature file check) and sets `varPtmTypeNum`
  from the modification file itself, and `TopPIC_topmg_identify` searches
  every spectrum with an empty `raw_prsm_ext`.
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

## ProForma strings (`seq/proteoform_util`, issue #610)

`proteoform_util::getProFormaStr` writes a proteoform in ProForma 2.1
(HUPO-PSI, https://github.com/HUPO-PSI/ProForma; the spec and its compliance
levels are in that repository's `SpecDocument/ProForma_v2.1_Final.pdf`).
`PrsmMatchTableWriter` puts it in the `ProForma` column after `Next amino
acid` (PrSM, proteoform and protein tables of toppic and topmg) and
`PrsmSqlWriter` in the `proforma` column of the SQLite tables; the manual
describes it in `docs/manual/proteoform.md` section 6. Facts to keep
straight (toppic_test experiment 128):

- **PTMs are Unimod name and accession**, `[Oxidation|UNIMOD:35]`,
  `[Cation:Na|UNIMOD:30]` (joined with `|` as synonyms, ProForma 2.1
  section 7.5, a Level 2 feature). The name is `PtmBase::getUnimodName(id)`,
  the abbreviation of the `ptm_base.xml` entry with that Unimod id, never a
  modification file's name (unchecked: `Sodium,...,30` gives
  `[Cation:Na|UNIMOD:30]`); an id without an entry (newer than the file)
  gives the accession alone. A PTM without a Unimod id (a modification
  file's line with UnimodID -1) is its mass with the name as an `INFO` tag
  (section 6.6), `[+12.3457|INFO:MyPTM]`. A name with `|`, `#` or
  unpaired brackets is left out (`isProFormaName`; `Cation:Mg[II]` is
  written). A label (`#s1(0.900)`, `#g1`) goes after the tag's first part
  (`labelTag`): after an `INFO` part it would be read as part of the
  information.
- **`ptm_base.xml` is Unimod**: `No PTM` (TopPIC's only own entry; the
  empty PTM is found by that name, not by mass 0) and every Unimod
  modification (1,560, unimod.obo of 17 Feb 2026) in id order, written by
  toppic_test's `tests/202610/scripts/proforma_output/gen_ptm_base.py`:
  abbreviation = the Unimod (PSI-MS) name, name = the Unimod full name,
  mass = the exact sum of the AME2020 nuclide masses over the Unimod
  composition (bricks expanded). Regenerate it for a new Unimod release;
  do not edit it by hand. Experiment 128 renamed the earlier entries to
  their Unimod names (Carbamidomethylation -> Carbamidomethyl, Oxoalanine
  -> Cys->Oxoalanine, Allysine -> Lys->Allysine, Thiocarboxy ->
  Carboxy->Thiocarboxy, Cholinephosphoryl -> Phosphocholine, TPQ ->
  Quinone; carboxylation and Citrullination were the same Unimod PTMs as
  Carboxy and Deamidated) and removed Decarboxy, which has no Unimod
  entry. The old names are not kept as aliases (the user's choice): a
  FASTA sequence written `C[Carbamidomethylation]` stops with "Unknown
  PTM", a modification file's old name defines a PTM of its own with the
  file's mass, and result XML written with an old name cannot be read
  (result files of other versions are refused anyway, see TopDiff). The
  C57 PTM is found by `getC57AbbrName()` ("Carbamidomethyl").
- **Every Unimod name is built in**, so a modification file's line whose
  name is a Unimod name takes the Unimod entry (`PtmBase::getPtmPtr`,
  documented in `docs/manual/proteoform.md` section 3); when the file's
  mass differs by more than 0.01 Da, a `WARNING:` line is printed once per
  name (on `std::cout`: toppic and topmg set the log level to errors, so
  `LOG_WARN` would not show).
- **Where each mass shift goes**, per alteration (`Alter`), on the
  alteration's residues, or its mass shift's when the alteration has none
  within them (TopMG writes an unexpected shift on one residue with an
  empty alteration range; TopPIC's alterations all have their shift's
  range): a protein
  N-terminal acetylation (PROTEIN_VARIABLE, residue 0, the Acetylation PTM)
  and the truncation N-terminal modification are N-terminal modifications
  `[...]-`; an alteration on one residue follows it (fixed, input, variable,
  or an unexpected shift on one residue); an alteration on several residues
  is a range `(...)[...]`, alterations with the same range sharing it (a
  variable shift's unlocalized PTMs); a local annotation (MIScore, `-B`)
  over several residues is written at its best site with `#s<n>(score)` on
  every site with a positive score (`LocalAnno::getScrVec` is indexed by
  the residue position in the proteoform when its length is the
  proteoform's); an empty or out-of-range alteration is an unknown
  position `[...]?`.
- **Ranges may not overlap or nest** in ProForma (TopPIC's `Proteoform`
  string nests them): the ranges are taken from the left, wider first, and
  one overlapping a taken range becomes a set of positions `#g<n>` with
  its first residue preferred.
- Delta masses: 4 decimals with a sign always (`+0.0000` too).
- **TopMG's variable shifts carry the aligned mass** (an `Oxidation` range
  of 16.2826 Da, a `Phospho;Acetyl;Acetyl;Acetyl` range of 206.0429 Da
  against 206.0980), and its `Proteoform mass` column includes it; the
  ProForma accessions stand for the PTMs' exact masses, so a reader's mass
  of a TopMG string is the theoretical one and differs from that column
  (by up to 0.42 Da on 480_F3_02). TopPIC's variable shifts have the PTM
  masses.
- Checked with pyteomics.proforma (5.0.1) on 480_F3_02 (toppic default
  setting, -B, an N-term variable PTM, topmg, PTMs in the FASTA): every
  string parses, has the proteoform's residues and, for toppic, the
  Proteoform mass within the rounding of the written masses (topmg: after
  the aligned masses above).

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

## Activation of MS/MS spectra (`ms/mzml/pw_ms_reader`, `res/base_data/activation_base.xml`)

With `-a FILE` (the default) `PwMsReader::parseActivation` takes each
MS/MS spectrum's activation from the activation cvParams of its first
precursor, through `getActivationName` (branch `ethcd_activation`, issue
#530, toppic_test experiment 132). Facts to keep straight:

- **The set of terms decides, not their order.** MS:1002631 (EThcD), or ETD
  (MS:1000598) with beam-type CID (MS:1000422, pwiz `MS_HCD`) or supplemental
  beam-type CID (MS:1002678) -> `ETHCD`; MS:1003182 (ETciD), or ETD with CID
  (MS:1000133) or supplemental CID (MS:1002679) -> `ETCID` (EThcD is tested
  first, so ETD with both kinds of CID is `ETHCD`); otherwise `ETD`, `HCD`,
  `CID`, `MPD` (MS:1000435) in that order; no known term -> the "Unable to
  read the activation method" error, as before. The reader used to keep the
  first of CID/HCD/ETD/MPD in file order, and the converters order the terms
  differently: msconvert writes CID, ETD, ECD, PQD, HCD, then the supplemental
  term (EThcD: ETD, beam-type CID, supplemental beam-type CID, read as ETD;
  ETciD: CID, ETD, supplemental CID, read as CID), ThermoRawFileParser the
  first reaction, then the supplemental one (EThcD: ETD, supplemental
  beam-type CID). Neither writes MS:1002631 or MS:1003182, so do not wait for
  those terms alone. The vendored pwiz (CV 4.1.182) knows all of them;
  updating it changes nothing here.
- **`ETHCD` and `ETCID` have the ion types of ETD** (C / Z_DOT in
  `activation_base.xml`): TopFD writes `ACTIVATION=ETHCD` or `ETCID` in the
  msalign file and C / Z_DOT in the SQLite table `ms2_spectrum`, and toppic
  and topmg read the name back (`ActivationBase::getActivationPtrByName`)
  and search EThcD and ETciD spectra as ETD ones. Searching b / y (or c / y
  for ETciD's CID) ions with c / z-dot would need several ion types per
  `Activation` (one N- and one C-terminal ion type now; used by
  `prsm/theo_peak_util`, the PRM spectra's `getN_BYShift`/`getC_BYShift`,
  the zero-shift fast search and the MCMC E-values).
- **An activation name is listed in several places**:
  `activation_base.xml`, the `-a` help text and validation of
  `topfd_argument.cpp`, `topdia_argument.cpp`, `toppic_argument.cpp` and
  `topmg_argument.cpp`, the `activationComboBox` items of the five GUI forms
  (`topfddialog.ui`, `topdiadialog.ui`, `toppicwindow.ui`, `topmgwindow.ui`,
  `topmergedialog.ui`; the item text is the option value), the `-a` rows of
  the four manuals and `res/README.md`. Change them together.
- **Results** (toppic_test experiment 132, amendment 1): on
  PXD019792_191018_Mmazei_R2_3_ETHCD_5sec (Orbitrap Fusion Lumos, msconvert
  3.0.24284) all 2,190 MS/MS spectra are ETHCD (main: ETD); the MS/MS
  masses, charges and headers and the feature files equal main's (the MS1
  msalign differs in 7 of 1,590 spectra, main against a second main run in
  9: EnvCNN and the greedy search are not exactly repeatable). toppic on
  main's TopFD output gives main's files; with `ACTIVATION=ETHCD`, or `-a
  ETHCD`, the same but the Fragmentation column (346 PrSMs, 206
  proteoforms, 66 proteins at 1% FDR against the M. mazei proteome). Test
  spectra in 12 activation cases (msconvert's and ThermoRawFileParser's
  EThcD and ETciD, the two terms alone, EThcD with HCD first, ETD, HCD,
  CID, CID with HCD, MPD) all get the expected activation and ion types;
  the other 15 mzML files of toppic_test's `data/` (34,143 MS/MS spectra,
  HCD or CID) keep theirs. Amendment 2 (`ETCID`): the three ETciD cases
  (msconvert's, ThermoRawFileParser's, the term alone) give ETCID and ETD
  with beam-type CID and CID gives ETHCD; toppic with the label ETCID (from
  the msalign or `-a ETCID`) writes the ETD run's files but the
  Fragmentation column. None of the local files has ETciD spectra.

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
`PeakTolerance::compPrecErrorTole(prec_mass)` checks in the
filters/searches/E-value code (`MsHeader::getFirstPrecErrorTolerance(ppo)`
until the Dalton option below), `refinePrecursorByMatchedIons`, the
localization shift checks
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

## Precursor error tolerance in Dalton (`para/peak_tolerance`, `--prec-dalton-error`)

Branch `prec_dalton_error` (toppic_test experiment 133, report 0133, issue
#563; merged into main on 2026-10-07). toppic's hidden `--prec-dalton-error <Da>` (argument
`precDaltonError`, empty by default; checked to be a positive number) sets
`PeakTolerance::setDaltonTolerance` on the precursor tolerance in
`PrsmPara`, and every precursor tolerance is then that value instead of
the precursor mass x the ppm of `-E`. The parameter block prints it in Da.
Facts to keep straight:

- **One rule for a precursor tolerance.** `compPrecErrorTole(prec_mass)` is
  prec_mass x ppo (no minimum) or the Dalton value; it serves the filters'
  and searches' precursor tests, `refinePrecursorByMatchedIons`, the
  E-values' candidate window and isotope offset (`CompPValueArray` keeps the
  tolerance as a member; its functions and `EValueProcessor`'s lost the ppo
  parameter) and the MCMC tolerance. `compStrictErrorTole(mass)` returns the
  Dalton value for every mass but 0: the zero PRM peak keeps the minimum
  tolerance (0.01 Da), which no ppm value ever changed.
  `MsHeader::getFirstPrecErrorTolerance(ppo)` and
  `Precursor::getErrorTolerance(ppo)` are gone, and
  `getFirstPrecMonoMassMinusWaterError` takes the tolerance in Da.
- **Unchanged:** the +-1 Da precursor offsets of the stages, `-e`, `-p` and
  the index. Without the option the 480 run is main's byte for byte (60
  files).
- **The refinement of the shift searches.** `refinePrecursorAndHeaderShift`
  steps over +-2 x tole and pairs fragments within tole + their own
  tolerance. With fewer than two pairs its median index is 0, so the first
  bin wins and the precursor moves by +2 x tole (main 03a0c37a: 264 of
  1,315 one-shift PrSMs of 480 at +2 x 10 ppm; +1.0 Da at 0.5 Da).
    - **Fixed on main by experiment 135** (e720aec2, issue #514): PsAlign and
      PsAlignMultiType call `refinePrecursorByMatchedIons` (the delta in
      [-T, T] that matches the most precursor-dependent ions within the
      fragment tolerance, 0 when none matches). On 480: 3,354 / 1,996 / 440,
      no shift PrSM at +2T. This branch's own fixes (amendment 1: bfa89f11,
      8a5a4582, 75c7ae13) were reverted (3413086c, amendment 2).
    - **TopMG's graph alignment** (topmg and toppic `--topmg`) kept
      `refinePrecursorAndHeaderShift` with the fallback until toppic_test
      experiment 143, which gave it `refinePrecursorByMatchedIons` too and
      removed the old function (section "Precursor refinement of TopMG's
      graph alignment").
    - **In the merge** `refinePrecursorByMatchedIons` takes the precursor
      tolerance instead of its ppo (`compPrecErrorTole`, the same T without
      the option), so `--prec-dalton-error d` sets T = d: the shift searches
      refine the precursor over +-d by the most ions matched with the
      fragment tolerance.
- **Shift E-values depend on the precursor peak width.**
  `CompProbValue::compPrecProbs` multiplies the one-shift precursor
  probability by the precursor peak's width and the i-shift one by the
  width again, so a wider tolerance lowers the E-values of PrSMs with
  unexpected shifts (median log10 change at 0.5 Da: -0.59 one shift, -1.16
  two shifts) and of the decoys among them; zero-shift E-values rise
  (+0.31), from the candidate window.
- **Results on 480** (setting d, 1% FDR; main 3,337 / 2,001 / 437):
  0.05 Da 3,334 / 1,997 / 439, 0.2 Da 3,322 / 1,985 / 443, 0.5 Da 3,316 /
  1,980 / 443; the E-value cutoff falls from 2.86 to 0.80 at 0.5 Da (shift
  decoys). 0.2 and 0.5 Da take 7% longer (zero-shift and variable PTM
  filters). With precursors moved by up to +-0.3 Da, 10 ppm reports 3,174
  PrSMs but only 295 equal to main's (the errors become unexpected
  shifts); 0.5 Da reports 3,057, of them 1,417 equal to main's (1,140 of its
  1,471 zero-shift PrSMs), with an E-value cutoff of 0.042.

## Graph alignment tolerance (`search/graphalign`, TopMG)

TopMG's graph alignment matches the mass differences (distances) of the
spectrum graph with those of the proteoform graph within
`GraphAlignMng::topmg_align_error_tolerance_` (Da; a constructor argument,
so both drivers must pass it), applied as `getIntTolerance()` in the
integer masses of `convert_ratio_` (274.335215 per Da: 0.1 Da is 28
units) by `graph_align.cpp` (distance matching and cutoffs),
`SpecGraph::compDistCover` and `graph_post_processor.cpp` (PTM masses this
close are merged). It is the argument `topmgAlignErrorTolerance`, default
0.1 Da (the former hard-coded `error_tolerance_`): topmg's hidden
`--align-error-tolerance` (visible in `-h` until 2026-10-08, when it was
moved to the hidden options and topmg_gui's field "Alignment error
tolerance (Da)" on the Advanced Parameters tab was removed, with its
`--align-error-tolerance` entry in `gui/util/command.cpp`, so the GUI leaves
topmg at 0.1 Da; the parameter block still prints it) and
toppic's hidden `--topmg-align-error-tolerance` for `--topmg pipeline` /
`topmg-only`, as `topmgShiftNumber` is `--num-shift` / `--topmg-num-shift`.
Both parsers refuse a value that is not a positive number. It is
independent of `-e` / `-E` (ppm): the precursor tolerance enters the graph
step only through the spectrum graph's zero and precursor PRM peaks
(`spec_graph_reader`) and the precursor refinement after an alignment.
The MCMC E-values (`MCMCMng::error_tolerance_`) keep their own hard-coded
0.1 Da. Facts to keep straight:

- **The post-processing keeps the variable PTM limit.**
  `GraphPostProcessor` explains each variable shift by a PTM combination of
  `mass_ptm_map_`, which holds combinations of up to `max_known_mods_` PTMs,
  a combination being left out when a mass within the tolerance is already
  there. Each shift takes its nearest combination; when these exceed
  `max_known_mods_` (`--var-ptm`) PTMs in all, `choosePtmCombs` chooses the
  combinations together: the smallest sum of mass differences with at most
  `max_known_mods_` PTMs (a dynamic programme over the shifts and the PTM
  number). Below 0.1 Da fewer combinations are left out, and the nearest
  ones gave proteoforms more variable PTMs than `--var-ptm` (6-9 of 5 at
  0.02 and 0.05 Da on the first 200 spectra of 480_F3_02; 4-5 of toppic's
  3), which overran `ptm_mass_vec2d` in `mcmc_dpr_processor.cpp` (indexed
  by the variable PTM number; a segmentation fault). Main 1f88b6e8 kept the
  limit by giving each shift in turn only the PTMs left, which labelled a
  70.06 Da shift of H2B1N as one Phospho (mass split 70.06); the joint
  choice gives Acetyl + 2 Methyl (toppic_test experiment 114). Wherever the
  nearest combinations keep the limit (every PrSM at 0.1 and 0.5 Da) they
  are the choice, so those outputs are unchanged.
- **TopMG's final tables are not repeatable with more than one thread.**
  `comp_pvalue_mcmc.cpp` shuffles with one global `std::mt19937 rng(0)`
  that all MCMC threads share, so the E-values differ between runs (1-2 of
  344 identical between two `-u 16` runs of main 0f938697, 110 and 113
  PrSMs); the intermediate files up to `topmg_graph_post` are repeatable
  and are what to compare (`-k`). toppic's `--topmg` stage with its default
  generating-function E-values is repeatable.

## Protein N-terminal forms in TopMG (`search/graph`, `search/graphalign`)

TopMG's graph alignment (topmg and toppic's `--topmg`, the same code)
reports the protein N-terminal forms of `-n` (`NONE`, `NME`,
`NME_ACETYLATION`, `M_ACETYLATION`) as toppic does: the "Protein N-terminal
form" column, the start residue, the `[Acetyl]-` prefix, the acetylation
not a variable PTM. Facts to keep straight:

- **Graph edges** (`ProteoAnno::anno`, only for the protein's N-terminal
  sequence and only the forms `prot_mod_util::allowMod` accepts for it):
  `M_ACETYLATION` is an acetylated-M edge at residue 0; `NME_ACETYLATION_X`
  is an edge at residue 0 that weighs the acetyl group alone (an empty
  amino acid with the Acetyl PTM: M removed, residue 1 acetylated). Both
  are `PROTEIN_VARIABLE` and take one of the alignment's `--var-ptm`
  layers. Being at residue 0, only an alignment that starts at the protein
  N-terminus can use them. The NME acetylation used to be an acetylated
  residue 1, open to an alignment that kept the M (`M(A)[Acetyl]...`).
  `NME` is the start at vertex 1; `is_nme_` (reset for every protein) says
  whether that start is the protein N-terminus (`getFirstDiagonal`), and
  with `-w` an alignment starts at vertex 1 only then.
- **The form is read back in `GraphPostProcessor`** (`chooseNTermForm`),
  before the variable shifts are named, from the start and the variable
  shift that begins at the first residue (the edge's layer change makes the
  first diagonal end at the start node, so that shift holds the edge's mass,
  possibly with PTMs of the next residues). Start 0: `NONE`, `M_ACETYLATION`
  (shift - acetyl) or `NME_ACETYLATION` (shift - acetyl + M, the
  proteoform then starts at residue 1: residues, fixed and other shifts
  moved by one); start 1: `NME` if allowed, else `NONE` (a truncation), or
  `NME_ACETYLATION` (shift - acetyl); later starts `NONE`. An acetylated
  form is taken when what it leaves (nothing, or the nearest PTM
  combination) explains the shift at least as well as the shift's own
  naming (error in Da within the tolerance, or no larger than the shift's
  own error up to twice the tolerance) with fewer variable PTMs (or as
  many and a smaller error). The form is built as
  `proteoform_factory::geneProtModProteoform` builds it: the `ProtMod`, the
  residue sequence's N-terminal mod (its mass), a `PROTEIN_VARIABLE` shift
  on residue 0 (the annotation). A mass left over (at most twice the
  tolerance) is dropped, so the proteoform mass can differ from the
  adjusted precursor mass by that much. Main used to name the acetylation
  as a variable `Acetyl` (adding Acetyl to the naming whenever an
  acetylated form was allowed; now only the `-i`/`-b` file's PTMs are
  names) and set `NME` for every start at residue 1.
- **Unexpected shifts are not read as forms**: an N-terminal acetylation
  that an alignment explains with an unexpected shift (topmg's default
  `--num-shift 1` allows one) stays one.
- **MCMC E-values** (`comp_pvalue_mcmc`): the residue sequence's N-terminal
  mod is added to the N-terminal fragment masses of the PrSM and of the
  random sequences (`n_term_shift_`); the candidate counts already have
  the forms (`CountTestNum`). TDGF (toppic `--topmg`) handles the forms as
  for toppic's own PrSMs.
- **Truncation N-terminal modifications** (proteoform N-terminal variable
  PTMs, Position `N-term` of the PTM file; toppic's since experiment 100,
  "Truncation N-terminal modifications" below). The graph has no edge for
  them: `GraphAlign::getConsistentPairs` also matches each distance of the
  spectrum's zero vertex (the N-terminus) with proteoform distances m
  shorter (m: each N-term mass, `GraphAlignMng::n_term_int_masses_`, filled
  by `GraphAlignProcessor::process`) from the vertices where the mod is
  allowed (`GraphAlign::n_term_starts_`: protein residue 3 or later, an
  amino acid of the line), one variable PTM layer more; the spectrum's
  distance cover marks those distances (`SpecGraph::compDistCover`).
  `chooseNTermForm` reads the mod back from the variable shift that begins
  at the first residue of a proteoform starting at residue 3 or later (the
  same rule as the protein forms) and sets it with
  `Proteoform::setTruncNTermModPtr`. Variable shifts are named only with the
  side-chain PTMs (`readModTxt` list 2), in the post-processing and in the
  MCMC's combinations (an N-term line used to name side-chain shifts, and
  one name on an N-term and a side-chain line made `readPtmTxt` throw); the
  MCMC candidate number of a PrSM with the mod sums over the N-term masses;
  `toppic --topmg topmg-only` sets `nTermVarPtmTypeNum` (the TopMG stage's
  and the cluster correction's generating function multiply by T_N + 1).
  **Ties go to the path without the mod**, so that a protein N-terminal
  form wins a tie against it: `GraphAlign` keeps a flag per DP entry
  (`trunc_paths_`; a mod's consistent pair stores its cell c as -1 - c),
  and wherever two paths score the same (the best predecessor, a variable
  against an unexpected-shift path, the best shift node, the end cell, and
  equal matched fragment numbers over the PTM numbers in `geneResult`) the
  one without the mod is kept. Without N-term lines nothing changes.
  Results (toppic_test experiment 116, report 0116; 480_F3_02, topmg
  `-f C57 -d -t FDR -T FDR`, oxidation M, phospho STY, `Acetyl,*,N-term`):
  topmg reports 12 and 16 truncation N-terminal acetylations at 1% FDR,
  8 and 9 of toppic's 11 (DDX17 from residue 25, JUPI2 44, CLCB 69, TRIR
  28), `--topmg topmg-only` 39 with 9 of them (main: none); PrSMs 1,973
  and 2,022 against 2,005 (MCMC noise), no time cost. The modification
  competes with the protein forms: the alignment keeps the path with the
  most matched fragments whatever its PTM number, so an acetylated
  truncation plus PTMs can displace a protein form (24 PrSMs of the
  ported TopMG start elsewhere than toppic's, 11 where toppic has NME or
  NME_ACETYLATION). With the tie rule (fcba5f51) topmg reports 14 (9 of
  toppic's) and `--topmg topmg-only` 26 (the same 9): where toppic has a
  protein form 6 instead of 11, the other unconfirmed ones 11 instead of
  19 (H2B2E, H2B2F, CLCB now toppic's forms with its E-values). The 6 left
  are not ties: the truncation matches at least as many fragments as the
  best path main finds, or (YBOX3, one spectrum) has the higher DP score in
  the protein form's layer, so the form is dropped before its fragments
  are counted; giving the mod's paths their own DP layers would let both
  reach the fragment count.
- **Results** (toppic_test experiment 115, report 0115; 480_F3_02, topmg
  `-f C57 -d -t FDR -T FDR`, oxidation M, acetyl K, phospho STY): on the
  spectra topmg and toppic both identify with the same protein and toppic
  starts at residue 1-2, the same form and start in 89.9% (main 21.7%;
  `--topmg topmg-only` 81.6%, main 18.5%); the final topmg table has 299
  NME_ACETYLATION and 35 M_ACETYLATION PrSMs and no variable Acetyl where
  an acetylated form fits (main: 188 such PrSMs, 140 Acetyls on ranges
  without K). topmg 2,009 and 1,995 PrSMs against 1,854 (the acetylated
  forms' MCMC E-values, a median 10^0.63 better); `--topmg topmg-only`
  2,976 against 2,979 (TDGF: one PTM fewer, x1/4). No time cost. The
  remaining disagreements with toppic are alignments: histone H2B starts,
  isobaric alternatives (acetyl-S N-terminus against a C-terminal E;
  acetylation + K adduct against Phospho), an acetyl put on a K further in.

## PTMs written in FASTA sequences in TopMG (`search/graph`, `search/graphalign`)

A FASTA sequence can carry PTMs in the extended format, a residue followed
by the PTM's abbreviation of `res/base_data/ptm_base.xml` in brackets
(`M[Oxidation]`, `C[Carbamidomethyl]`; the names of earlier versions, such
as `C[Carbamidomethylation]`, are no longer read, see "ProForma strings").
`FastaSeq::compAcidPtmPairVec` takes a name up to the `]` that closes its
`[`, counting nested pairs, so the 51 Unimod names with brackets can be
written (`E[Cation:Mg[II]]`; it used to stop at the first `]`); an
unpaired `[` stops with "Unpaired '[' after residue <n> in the sequence
of <protein>".
`FastaSeq` reads them into residue/PTM pairs (`getAcidPtmPairVec`), and
`proteoform_factory::geneDbProteoformPtr` builds the modified residues with
an `INPUT` mass shift on each (shown in the Proteoform column as
`(M)[Oxidation]` and in the Database protein sequence; not in the Fixed
PTMs column or the variable PTM counts). Since branch `topmg_fasta_ptm` TopMG's graph alignment (topmg and
toppic `--topmg`) treats them as fixed modifications. Facts to keep
straight:

- **The mass graph is built from the pairs** (`ProteoAnno::anno` takes a
  `StringPairVec`). It used to take the raw string, whose
  `convertStrToResiduePtrVec` dropped the brackets and lower-case letters
  and kept the rest: `M[Oxidation]` became M and O (pyrrolysine), the graph
  had more vertices than the proteoform had residues, and topmg and toppic
  `--topmg` ended with a segmentation fault at the start of the graph
  alignment (toppic_test experiment 122).
- **One edge per modified residue.** A residue with a PTM, of the FASTA or
  of `-f`, is marked `FIXED`, and `getMassGraphPtr` gives a position with a
  `FIXED` residue that edge alone (and `PROTEIN_VARIABLE` edges, below): no
  variable PTM edge. A first residue with a FASTA PTM gets no
  `M_ACETYLATION` edge, which would modify it, and
  `GraphPostProcessor::chooseNTermForm` does not take that form for it. The
  forms that remove it stay: NME (a start, not an edge) and
  `NME_ACETYLATION`, whose edge in place of the M is the position's
  `PROTEIN_VARIABLE` edge; the post-processing drops the FASTA PTM of the
  removed M (an `INPUT` or `FIXED` shift on residues before the form's
  start). A first version took every acetylated form away from such a
  residue: with M1 written oxidized, YBOX3, FUBP2 and H1.5 lost their
  NME_ACETYLATION PrSMs on 480_F3_02 (experiment 122). The truncation
  N-terminal modifications (not edges either) stay open, as for a fixed
  modification's residue. A `-f` modification is not applied to a residue
  with a FASTA PTM (`applyFixedMod` compares with the unmodified residue).
- **The post-processing keeps the `INPUT` shifts** (it rebuilds each
  proteoform from its `FIXED`, `VARIABLE` and `UNEXPECTED` shifts, which
  dropped them from the reported proteoforms).
- **Residues are registered before the threads start.**
  `fasta_util::dbPreprocess` ends with `addPtmResidues`, which adds the
  residues of the database's FASTA PTMs to `ResidueBase` (a bracket is
  looked for first; a database without PTMs is read once more as text and
  nothing else happens). `ResidueBase::getBaseResiduePtr` appends a new
  residue to a shared vector without a lock, which the search threads did
  for a PTM not in the `-i`/`-b` file. An unknown PTM name now throws there.
- **E-values**: `count_util::updateResidueCounts` (the residue frequencies
  of `CountTestNum`: TDGF and the MCMC's candidate numbers) adds a residue
  that is not in the list, so that a residue with a FASTA PTM is counted as
  its own residue, as a fixed modification's residue is; it used to be
  left out (toppic `--topmg` E-values a median 10^0.1 larger with every C
  written as `C[Carbamidomethylation]` than with `-f C57`). The MCMC's
  random residues (`mcmc_mass_table_util`) take the fixed modifications of
  `-f` only.
- **Protein hits** (`prsm/search_fasta_match`, the other proteins of a PrSM
  in the TSV and SQLite tables of toppic and topmg) index the residues, one
  letter each, and a hit must have the same FASTA PTMs. The raw string was
  indexed and cut at the proteoform's residue positions: with brackets
  before the start the cut was a piece of a PTM name and the next residues,
  which matched thousands of proteins (350,183 rows of other proteins in
  topmg's PrSM table of 480_F3_02 with every C written as
  `C[Carbamidomethylation]`, against 269 with `-f C57`). A database
  without PTMs is indexed and searched as before.
- The PTM localization of toppic (`-l`, `stat/local/local_proteoform`)
  still reads the raw string (`getRawSeq`); topmg has no such stage.
- Results (toppic_test experiment 122, report 0122; 480_F3_02, `-d -t FDR
  -T FDR`, oxidation M, phospho STY, `Acetyl,*,N-term`): every C written
  as `C[Carbamidomethylation]` without `-f` gives `-f C57`'s filter,
  alignment and post-processing files (topmg) and all of `--topmg`'s files
  with its E-values and tables (but the Fixed PTMs and Database protein
  sequence columns), with the same index (built in memory; the shared
  index files of toppic_test's `data/`, built on 2026-09-28, give 70 of
  60,135 one-PTM filter candidates differently). Without FASTA PTMs the
  files are main's. With the 249 methionines main reports oxidized written
  as `M[Oxidation]`, 86-90% of the PrSMs that had the oxidation are the
  same proteoform with it from the FASTA (`--topmg` E-values a median
  10^-0.48 smaller), and the unoxidized forms over those sites take an
  unexpected shift or are lost (`--topmg` 2,841 PrSMs against 2,959).

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

## E-values of PrSMs found at a +-1 isotope offset (`stat/tdgf/comp_pvalue_array`)

Branch `evalue_adjusted_prec` (toppic_test experiment 117, report 0117;
merged into main on 2026-10-05). `CompPValueArray::process` (the shared mode that every toppic
E-value pass uses: the first pass, the multiple-shift pass, the cluster
protein correction and `--topmg`'s TDGF E-values) groups the PrSMs of a
spectrum by isotope offset. Facts to keep straight:

- **The offset rule** (`getIsotopeOffset`) is the one the candidate number
  already used. A PrSM without unexpected shifts, variable PTMs or a
  truncation N-terminal modification has offset +1 / -1 when its adjusted
  precursor mass is more than the precursor error tolerance above / below
  the spectrum's label-free precursor mass. All other PrSMs have offset 0,
  including the one-shift and variable-PTM offsets of the opt-in
  `--enable-*-prec-error` options, as on main.
- **Offset 0** keeps main's computation: the spectrum set's PRM spectra
  and generating function. Their E-values are unchanged (equal on every
  PrSM of 480_F3_02); a smaller maximum score in the group only changes the
  height of the generating function's table, whose top row saturates.
- **Offset +-1** gets `prm_ms_factory::geneMsSixPtrVec` for the observed
  precursor mass +- 1.00228 Da, the PRM spectra the zero-shift search built
  for that mass, so the generating function's precursor mass is the
  adjusted one. `compMultiExpectedValues` takes the candidate number's
  offset from the PRM spectra's precursor mass instead of the PrSM's
  original precursor, so that the offset is not applied twice. This also
  fixes the separate mode (`compSingleExpectedValue`, unused by toppic),
  which already built the PRM spectra at the adjusted mass and then shifted
  the candidate number once more.
- **Results on 480_F3_02** (`-s 2`, oxidation only):
    - For 94 spectra that experiment 116's oracle searched with the right
      precursor, the new offset E-values equal the direct ones exactly. On
      main the same matches differ by 10^-16.9 to 10^2.6.
    - 1% FDR: 3,336 / 1,996 / 437 against main's 3,337 / 1,997 / 438 on
      TopFD 47b9fc57's output, 3,339 / 2,000 / 437 against 3,340 / 2,001 /
      438 on TopFD ead67ddc's. The E-value pass takes 4-5 s longer (+4%).
    - Decoy E-values move little (log10 ratio within -1.35..1.12); target
      offset PrSMs move by up to 17 orders of magnitude, in both directions.
    - **The best PrSM changes for deamidated spectra.** For an Asn-Gly
      deamidated proteoform whose precursor is right, the unmodified
      zero-shift match through the -1 offset now gets the E-value a right
      precursor gives and wins over the deamidation PrSM (1 shift, (N)[+0.98]),
      with 4-11 fewer matched fragments. Deamidated PrSMs at 1% FDR fall
      from 50 to 45 and from 60 to 53.
    - **Deamidation as a variable PTM repairs only part of it**
      (`Deamidated,0.984016,NQ,any,7` next to oxidation, report 0117
      section 8): 86 / 96 PrSMs with a Deamidated variable PTM against
      main's 91 / 102. In HMGN1 scans 2833, 2838 and 2974 the offset match
      still wins: the generating function on the PRM spectrum of P - 1
      gives its fewer fragments a smaller probability per candidate than the
      one on P gives the variable-PTM PrSM, although the precursor mass is
      4-7 mDa above the deamidated mass and 11-15 mDa below the unmodified
      mass + 1 isotope. Not implemented: choosing by precursor error between
      a PrSM at P and an offset PrSM of the same protein and positions.

## E-values of the variable PTM PrSMs with modified residues (`stat/tdgf/comp_pvalue_array`, `--disable-var-ptm-evalue-mod-residues`)

Branch `var_ptm_evalue_mod_residues` (toppic_test experiment 129, issue
#334), merged into main on 2026-10-06 with the method on by default. The
PrSMs of the variable PTM search get their own generating function; the
hidden toppic option `--disable-var-ptm-evalue-mod-residues` (argument key
`varPtmEvalueModResidues`, default `"true"`; the parameter block prints
"Variable PTM E-values with modified residues: false" only when it is set)
gives them the residues without variable PTMs, as before the merge. Facts
to keep straight:

- **The residues.** Without the method (the other PrSMs, and every PrSM
  before the merge) the generating function draws random sequences from
  the 20 residues with the fixed modifications (C57), each with its
  frequency in the search database (`CountTestNum`), and, for COMPLETE and
  PREFIX PrSMs, the first residues of the protein N-terminal forms. With
  the method, a residue that k - 1 variable PTMs can modify becomes k
  residues, itself and its modified forms, which share its frequency x,
  x / k each (`count_util::geneModResFreqs`): with oxidation on M, M ->
  M and M[Oxidation] at x / 2; with oxidation M, acetyl K, methyl KR and
  phospho STY, K -> K, K[Acetyl], K[Methyl] at x / 3. The PTMs are the
  "any" lines of the `-b` file (`mod_util::readAnywhereModTxt`, set in
  `TdgfMng::var_ptm_mod_ptr_vec_` by `setVarPtmModResidues` in
  `toppic_process.cpp`); N-term lines are not used.
- **Which residue a PTM modifies** follows the variable PTM search
  (`VarPtmAlign::initAllowShiftTable`): the residue must equal the
  modification's original residue (`Residue::isSame`, the amino acid
  without a PTM). A C57 residue takes no variable PTM. The N-terminal
  residue list is split the same way; its residues carry no protein
  N-terminal acetylation (that is an alteration of the proteoform, not of
  the residue), so an acetylated first M splits too, as the search lets it
  be oxidized.
- **Which PrSMs.** `CompPValueArray::useModResidues`: no unexpected shift,
  and variable PTMs or a truncation N-terminal modification, the PrSMs
  whose candidate number already takes the residue sequence mass.
  `process` splits the offset-0 PrSMs of a spectrum into these and the
  others (the +-1 offset PrSMs never are variable PTM PrSMs) and computes
  the generating function once for each group, so a spectrum with both
  computes it twice; the variable PTM group has no unexpected shift and
  needs only layer 0. The other PrSMs' E-values are unchanged: their group
  has a smaller or equal maximum score, which only lowers the saturated
  table height (as in the offset groups above).
- **One table.** `CompProbValue::setResidues` replaces the residue masses,
  frequencies and average residue length before a computation, so one
  CompProbValue (138 MB page table per thread) serves both lists. The
  average residue length of the new list (`count_util::computeAvgLength`)
  is used with it.
- **The candidate number is unchanged**, including the factor
  (variable PTM types + 1) per variable PTM.
- **Passes.** toppic's first E-value pass, the multiple-shift pass and the
  cluster protein correction get the PTMs; TopMG hosted in toppic
  (`--topmg`) does not (its PrSMs come from the graph alignment), and
  neither does topmg.
- **With `--disable-var-ptm-evalue-mod-residues`** toppic gives the files of
  main before the merge (4ecb19ea) byte for byte (56 of 56 files of a 480
  run).
- **Results on 480_F3_02** (`-s 2`, TopFD of main 4ecb19ea, toppic_test
  report 0129):
    - Oxidation only (`var_mods.txt`): 3,337 / 2,001 / 437 against
      4ecb19ea's 3,342 / 2,000 / 437, one decoy oxidation PrSM crossing the cutoff. The
      variable PTM PrSMs' E-values move little (median -0.004 log10 at 5-9
      matched fragments, up to +2.5 for ACBP with 48 fragments).
    - Four PTMs (oxidation M, acetyl K, methyl KR, phospho STY): 3,377 /
      2,013 / 442 against 3,371 / 2,012 / 442, but variable PTM PrSMs 731 ->
      800 and zero-shift PrSMs 1,469 -> 1,425. Their E-values fall (median
      -0.27 log10 at 5-9 fragments, -1.2 at 30 or more, decoys alike), because
      splitting S, T, Y, K, R and M moves 16% of the residue probability away
      from the unmodified residues. The PrSMs with variable PTMs then win ties
      and near-ties against the zero-shift PrSMs of the same spectrum, which
      keep the residues without PTMs. In 19 scans a protein N-terminal acetylation
      becomes an acetylated lysine near the N-terminus, with fewer (12) or
      equal (7) matched fragments. One H2B cluster flips from H2B2F to H2B1C
      + methyl (same mass; the cluster protein correction relabels 18
      scans).
    - The first E-value pass takes +6% (oxidation) and +18% (four PTMs), the
      whole run +1.0% and +1.6%; memory is unchanged.

## Precursor refinement of the variable PTM search (`search/diag/diagonal_util`, `search/varptmsearch/var_ptm_align`)

Branch `var_ptm_refine_prec` (toppic_test experiment 121, report 0121;
merged into main on 2026-10-05). Facts to keep straight:

- **The headers.** `VarPtmSlowMatch::geneVarPtmNTermShiftHeaders` builds one
  diagonal header per variable PTM mass combination, with C-terminal shift =
  observed (label-free) precursor mass - sequence mass - N-terminal shift.
  The C-terminal diagonal is the one whose C-terminal shift is within the
  precursor tolerance, so its C-terminal shift is the precursor error.
- **The refinement.** `VarPtmAlign::geneResult` gets `refine_prec_mass`
  (the proteoform mass plus the variable PTMs of the C-terminal diagonal)
  from `diagonal_util::refineVarPtmPrecursorAndHeaderShift`. That function
  also sets every backtracked header's C-terminal shift from it
  (`DiagHeader::initHeader`, as the slow match built them).
  `diagonal_util::refineVarPtmHeadersBgnEnd` then generates C-terminal ions
  without the precursor error.
- **On main** the headers keep the observed precursor's C-terminal shifts,
  so every y ion of the range refinement carries the precursor error.
    - Scan 5211 of 480_F3_02 is a deamidated proteoform whose precursor
      TopFD's fragment pairs moved by an isotope, 0.0202 Da above the
      proteoform mass. Every y ion below 2,020 Da was lost.
    - The variable PTM shift range spread over the lost ions (report 0120),
      and no ion inside a shift's range is counted.
- **The shift searches** (`PSAlign`, `PSAlignMultiType`, TopMG's graph
  alignment) do their own precursor refinement
  (`refinePrecursorByMatchedIons`, the most precursor-dependent ions within
  the fragment tolerance; until experiments 135 and 143
  `refinePrecursorAndHeaderShift`, the median fragment error). The
  variable PTM search needs none: its proteoform mass is exact once the
  alignment has chosen the variable PTMs.
- **Unchanged:** the dynamic programming (`VarPtmAlign::dp`) and the
  first-pass PrSMs without variable PTMs (31,702 identical on 480_F3_02).
- **Results on 480_F3_02** (`-s 2`, TopFD ead67ddc output):
    - Oxidation: 3,336 / 1,998 / 437 against main's 3,333 / 1,998 / 437.
      Oxidation + deamidation: 3,336 / 1,992 / 436 against 3,331 / 1,989 /
      436. Decoys at fixed E-values are unchanged.
    - First-pass variable PTM PrSMs of the same spectrum and proteoform: 708
      / 743 gain fragments, none lose; E-values better in 250 / 383, worse in
      none. Most gains are at precursor errors of 0.015 Da or more.
    - Reported PrSMs: 15 / 16 of the same proteoform gain fragments (ACBP
      scans 2000-2008 +9 to +13; 5211 38 -> 52, E 9.2e-17 -> 4.2e-33).
      1 / 4 change proteoform (a one-shift PrSM becomes a variable PTM
      PrSM). 3 / 5 PrSMs are new, all borderline (E 0.006-2). None is lost.
- **Remaining (fixed by the next section):** the range of a variable PTM
  shift was still the gap between the refined diagonals, wider than the
  modifiable residues (5211's oxidation of M44 on P43-M44).

## Precursor refinement of the shift searches (`search/diag/diagonal_util`, issue #514)

Branch `c_term_prec_adjust` (toppic_test experiment 135, report 0135; not
merged). Facts to keep straight:

- **Why the refinement matters.** A C-terminal ion of mass m enters the
  alignment as the PRM P - (m - C shift), P the observed precursor. In a
  diagonal anchored at the N-terminus (`n_strict && !c_strict`) it is matched
  with ppo x (m + P), the fragment tolerance plus the precursor's error
  (`prm_ms_factory`, `compDiagPair`); a C-anchored diagonal does the same to
  its N-terminal ions. The PrSM (`Prsm::initScores`) matches the proteoform's
  ions with ppo x m. The refined precursor decides whether the ions the
  alignment matched with the relaxed tolerance (the precursor-dependent
  ions: the C-terminal ions of every diagonal but the last and the
  N-terminal ions of the last) still match.
- **`refinePrecursorByMatchedIons`** (PsAlign, PsAlignMultiType) moves the
  precursor by the delta in [-T, T] (T = P x the precursor ppo, or the Dalton
  value of `--prec-dalton-error`) that matches
  the most precursor-dependent ions within the fragment tolerance (each
  theoretical ion once over the spectra of the group). Within the best
  delta interval it takes the median error of its ions, clipped to the
  interval; of equal intervals the one closest to 0; delta 0 when no ion is
  in reach. The shifts move as before: the last diagonal's N-terminal shift
  and the other diagonals' C-terminal shifts.
- **`refinePrecursorAndHeaderShift`**, main's method until then: the median
  of 0.005 Da bins over +-2T, which can leave the precursor tolerance, and
  +floor(2T / 0.005) x 0.005 Da with fewer than two pairs (report 0133).
  It stayed for TopMG's graph alignment until toppic_test experiment 143
  and is removed (next section).
- **The variable PTM search needs neither:** its proteoform mass is exact
  (`refineVarPtmPrecursorAndHeaderShift`).
- **Results on 480_F3_02** (`-s 2`, TopFD main 38ba67dc output):
    - First-pass target one-PTM PrSMs: 26,805 C-terminal ions matched with the
      relaxed tolerance, 81.1% of them within the fragment tolerance at P;
      the PrSM keeps 92.8% with main's refinement and 94.9% with this one
      (multiple-PTM: 86.1% and 89.3%). Matched fragments on the same
      alignments +0.9% (decoys +2.4%).
    - Main's fallback hit 251 of 1,281 reported one-shift PrSMs (+2T, about
      +20 ppm in the precursor and the shift mass); none now.
    - With every precursor moved by up to +-5 ppm, 80.1% of 6,779 target
      one-PTM PrSMs get exactly the unmoved run's adjusted precursor (main's
      bins: 0.33 ppm median).
    - 1% FDR: 3,354 PrSMs, 2,001 proteoforms, 438 proteins (main 3,337 / 2,001
      / 437); E-value cutoff 2.49 (2.86); decoys at E <= 1: 26 (19). Same
      running time.
- **Unchanged:** the alignment (it still counts the ions only the relaxed
  window matches), the E-value's generating function (relaxed tolerance on
  the observed precursor).

## Precursor refinement of TopMG's graph alignment (`search/graphalign/graph_align`, `search/diag/diagonal_util`)

Branch `topmg_prec_matched_ions` (toppic_test experiment 143, report 0143;
merged into main as 5c9e2e3d on 2026-10-07, after `filter_pairs`; on the first
500 spectra of 480_F3_02 toppic in setting d writes the same files as main
20d5895b, apart from the tables' start and end times). `GraphAlign::geneResult(s, m)`, shared by topmg and toppic `--topmg`,
calls `refinePrecursorByMatchedIons` like `PsAlign` and `PsAlignMultiType`;
`refinePrecursorAndHeaderShift`, its two helpers and
`GraphAlignMng::refine_prec_step_width_` are removed. Facts to keep straight:

- **Where TopMG's precursor error enters.** The spectrum graph's vertices are
  PRM peaks of the observed precursor P (C-terminal ones P - (m - C shift)),
  matched by distance within the alignment error tolerance (0.1 Da); every
  diagonal's C-terminal shift is set from P (`geneHeaders`). The first
  diagonal's shift comes from the zero vertex, the last one's from the
  precursor vertex, an internal one's is the mean of its nodes' shifts.
- **Single diagonals.** An alignment without unexpected shifts and variable
  PTMs is one diagonal, anchored at the zero vertex: its N-terminal ions do not
  move with P. `getDepIonPairs` takes its C-terminal ions as the dependent ions
  and the refinement moves its C-terminal shift (be5fb6ff). With two or more
  diagonals the rule is unchanged (the C-terminal ions of every diagonal but
  the last, the N-terminal ions of the last). toppic's searches never pass one
  diagonal: the one-PTM search asks `PsAlign` for one shift, the multiple-PTM
  search for two or more.
- **Downstream.** The PrSM's adjusted precursor is the refined one;
  `GraphPostProcessor` and the MCMC take it, the generating function of toppic
  `--topmg` works on the observed precursor's PRM spectra (its E-values did not
  change when only the single diagonals' precursor moved).
- **Results on 480_F3_02** (topfd main 7b5112c2, `-f C57 -d -t FDR -T FDR`,
  Oxidation M, Phospho STY, Acetyl N-term; toppic `--topmg topmg-only`):
    - The old method put 37.6% of 9,885 target graph alignments at its +2T
      fallback (decoys 50.9%; 60% of the alignments with fewer than 5
      fragments, 51% of single diagonals) and 44.3% outside [-T, T]; 281 of
      2,981 reported PrSMs at +2T. Now none outside T.
    - Same alignments: 1,024 of 9,238 targets gain matched fragments, 170 lose.
    - Reported PrSMs without unexpected shifts: adjusted precursor within 1 ppm
      of the proteoform mass for 96.1% (old 54.9%), > 10 ppm off 0.2% (6.0%).
    - Every precursor moved by up to +-5 ppm: 82.1% of the reported PrSMs get
      exactly the unmoved run's adjusted precursor (old 0%); strong alignments
      with one diagonal 98.9%, two 90.5%, three or more 35.0%.
    - 1% FDR: 2,968 / 1,748 / 403 PrSMs / proteoforms / proteins against
      2,981 / 1,754 / 398; more targets at every fixed E-value (E <= 0.01:
      2,752 against 2,736) and 4 / 10 more decoys at E <= 10 / 100, so the
      cutoff fell from 55.3 to 24. On the moved precursors 2,957 / 1,744 / 405
      against 2,954 / 1,738 / 402. Same running time.
    - topmg (MCMC E-values, two runs each): 1,979 / 1,965 PrSMs against
      2,021 / 2,002, proteoforms 1,212 / 1,217 against 1,239 / 1,219,
      proteins 390 / 390 against 394 / 386. The refined decoys gain fragments
      as the targets do (the old fallback hit half of the decoy alignments):
      decoys at E <= 0.1 114 / 121 against 93 / 98. Unchanged PrSMs' MCMC
      E-values scatter as between two runs of one binary (sd 0.9 in log10).
      MCMC +1.6 min (15,423 post-processed PrSMs against 14,818).
- **Known effects:**
    - Three or more diagonals recover the precursor for only 35%: an internal
      diagonal's mean shift mixes N-terminal and C-terminal PRM vertices, so it
      carries part of P's error.
    - A variable-PTM diagonal that matched no ion was dropped by
      `refineHeadersBgnEnd` and its mass merged into the unexpected shift; when
      the refined precursor lets one ion match, the diagonal stays and the
      unexpected shift can exceed `max_ptm_mass_` (500 Da), so the PrSM is
      rejected (spectrum 712 of 480, THOC4, E 2.5e-8 with the old method).

## The filters' mass list (`ms/factory/prm_ms_util`, `getIntMassErrorList`)

Branch `filter_pairs` (toppic_test experiments 138 and 142, reports 0138 and
0142; merged into main on 2026-10-07 as 5fb46b2c). Facts to keep straight:

- **What the list is.** The PRM peaks of a spectrum (prefix list: ORIGINAL
  m - N shift, REVERSED P - (m - C shift); suffix list: the same minus H2O)
  as integer masses round(x x `filter_scale_`, 100 = 0.01 Da) with the
  tolerance ceil(tol x 100) (ORIGINAL ppo x m, REVERSED ppo x (m + P)), plus
  the zero and precursor entries, sorted. The users add 1 per entry to every
  column in its window: `MassMatch::compScores` (one-shift and variable PTM
  filters), `compDiagMaxScoresImpl` (multi-shift filter) and
  `CompShiftLowMem::findBestShift` (the diagonal finder of `PtmSlowMatch`).
  The zero-shift filter does not use it; TopMG uses it in its one-shift
  filter (and the multi-shift one with `--use-asf-diagonal`).
- **One entry per integer mass and base type.** Peaks of one type with the
  same integer mass give one entry (the largest tolerance); an ORIGINAL and a
  REVERSED peak of the same integer mass are two. Main kept one entry per
  integer mass: at an accurate precursor a b ion and the PRM P - y of its
  complementary y ion differ by (precursor error - fragment error) x the
  masses, a few mDa, and often shared an integer mass, so the right
  proteoform scored one point for the pair (two when a bin boundary fell
  between the peaks). On 480_F3_02 an identified spectrum lost 10.9 entries.
- **Repeated integer masses are fine for the users**: they loop over the
  entries, and `findBestShift` takes its range from the first and the last
  entry (the precursor entry). Main already had them for combined spectra
  (the merging restarts per spectrum).
- **Results** (toppic_test experiment 142, against main 7b5112c2):
  480_F3_02 setting d 3,369 / 2,003 / 441 against 3,354 / 1,996 / 440, more
  targets at every rank of the target-decoy list; the TMT E. coli set of
  experiment 80 444 / 263 / 69 against 454 / 244 / 67 (more targets at most
  ranks; a lower 1% cutoff); TopMG in toppic (`--topmg topmg-only`) 3,025 /
  1,768 / 422 against 3,010 / 1,765 / 420. No change where the precursor and
  fragment errors differ by 10 ppm (P1_high_high_85min_1: nothing merges) or
  on st_1. No time or memory cost. Right proteoforms near a filter's top-N
  cutoff can still drop out when other candidates gain (4 of 23 changed
  spectra on 480).
- **The rule is a hidden option**, `--filter-mass-list
  <merge|fragment|single>` of toppic and topmg (argument key
  `filterMassList`, default `fragment`; branch `filter_mass_list`,
  toppic_test experiment 142, amendment 2; merged into main on 2026-10-08
  as abe33295, after `second_align`, `ps_align_penalty` and `local_align`).
  `SpPara::FilterMassListMode`
  (`kMerge`, `kFragment`, `kSingle`; `setFilterMassListMode` from the
  name, `isFilterMassListModeName` for the parsers' check) is set by
  `PrsmPara` and passed to `getIntMassErrorList` by its four callers, so a
  new caller must pass `getFilterMassListMode()` too. `merge` is main's
  rule up to 7b5112c2 (all peaks in one slot), `fragment` the rule above,
  `single` one entry per peak. Both parsers refuse other values, and the
  parameter block prints `Filter mass list:` only for `merge` and
  `single`, so default output is unchanged. toppic's `--topmg` stages take
  it from toppic's argument map. Checked on st_1 and 480_F3_02 (setting d;
  every output file but the time lines and that parameter line, 64 per
  run): the default and `fragment` write main c7f13d66's and filter_pairs
  57e44368's files, `merge` main 7b5112c2's, `single` filter_all_masses
  f6572a10's; topmg's one-shift filter files on st_1 equal those of the
  three binaries mode by mode.
- **`single` is not better** (toppic_test experiment 142, amendment 1,
  against `fragment`): its extra entries (61 per identified 480_F3_02
  spectrum, +13.6%) are 99% one fragment mass at its other charge states,
  which TopFD lists separately, so a fragment counts once per charge state.
  480_F3_02 3,357 / 2,000 / 437 against 3,369 / 2,003 / 441 (fewer targets
  at most ranks, strong PrSMs lost at the filters' cutoffs); TopMG in
  toppic fewer targets at every rank; topmg 2,017 against 2,052; TMT 525
  against 444 PrSMs, mostly because two high-ranked decoys are not found;
  P1_high_high 708 against 689 (a looser cutoff). The multi-shift filter,
  whose time grows with the square of the list, takes 34% longer on 480 and
  TMT, the runs 5-8%.

## Tag filter (`filter/tag`, `--tag-filter`)

Branch `tag_filter` (toppic_test experiment 150, report 0150; from main
d56ebc44, not merged). A sequence-tag protein filter whose candidates
(simple PrSMs) join the candidate lists of toppic's identification stages.
Off by default; the hidden toppic option `--tag-filter <stages>` (argument
`tagFilterStages`, a comma-separated subset of `var`, `one`, `multi`,
`topmg`, or `all`; any other name is refused) turns it on for the named
stages, with the hidden `--tag-filter-length <n>` (`tagFilterLength`, the
residues of a tag, default 4, at least `TagIndex::K`) and
`--tag-filter-cand-num <n>` (`tagFilterCandNum`, the proteins kept per
spectrum set, default 10). The parameter block prints the three only when
the filter is on. Facts to keep straight:

- **Tags** (`TagFilter::geneTags`). The PRM peaks of a spectrum set
  (`prm_ms_util::getPrmPeakPtrs` on `getMsTwoPtrVec`: each fragment mass
  as the N-terminal PRM, ORIGINAL, and as the precursor minus the mass,
  REVERSED, plus the zero and precursor peaks), per base type, peaks within
  the sum of their strict tolerances of each other merged (one fragment at
  several charge states). Two peaks whose difference is a residue mass (the
  19 letters with the fixed modifications; I is read as L) within the sum of
  their strict tolerances are an edge, one per residue that fits (K and Q
  both at a wide tolerance). A tag is a path of `tag_len_` edges in
  increasing mass, so it reads the protein N to C whichever series it comes
  from; the strings are looked up once each. About 180 tags per 480_F3_02
  spectrum (1,100 without the merging of charge states).
- **The index** (`TagIndex`). The search database (targets and decoys, the
  file the searches read) as one text with I as L, the positions of every
  4-mer in a bucket per 4-mer (as `prsm/search_fasta_match`), and the prefix
  mass before every residue (`residue_util::convertStrToResiduePtrVec` with
  the fixed modifications, so FASTA PTMs count; one double per residue,
  180 MB for the human database). A tag is looked up by its rarest 4-mer
  and verified at each position; every occurrence in every protein is
  returned. Built once per stage, 2-3 s.
- **Score** (`TagFilter::scoreHit`). An occurrence lies on the diagonal
  s = (PRM mass of the tag's first peak) - (prefix mass before the tag);
  the occurrences of a protein sorted by s form a diagonal while
  neighbours are within the sum of their first peaks' N-terminal
  tolerances (`getNStrictCRelaxTolerance`: the precursor's error enters a
  REVERSED peak's PRM). A diagonal's support is the number of distinct
  peaks its tags cover; the score is the best diagonal's support plus the
  second's when that reaches `TagFilterMng::getMinScore()` (tag_len + 3,
  default 7: a ladder of six residues or two tags apart on one diagonal),
  and a protein below the minimum is not a candidate. The first version
  counted the peaks of all occurrences: titin, nesprin and the other
  longest proteins collected hundreds of chance 4-mers and filled every
  list (37% decoys in the var lists, 7 of 174 confident zero-shift
  proteins of the first 500 spectra listed); with the diagonals and the
  minimum 83 of 174 are listed (one stage 88), 0.4% / 33% decoys. Ties go
  to the shorter protein, then the name (a name tie put every `DECOY_`
  before `sp|`). The top `cand_num_` proteins are the candidates.
- **What a candidate carries, by stage** (`TagFilterMng::Stage`; the
  spectrum sets of the stage's precursor error list, `getPrecErrorVec`, so
  the search pairs them by precursor mass as it pairs the filter's):
    - `one`: the diagonals' shifts (the best supported first, at most
      `max_shift_num_` = 10, two within 0.01 Da once) as `n_term_shifts`
      and, as c = P - protein mass - s (`getNTermShiftsFromCTermShifts`
      inverts it), as `c_term_shifts`; the same candidate in the four type
      files (`OnePtmSlowMatch` reads the lists its type needs).
    - `var`: `geneVarPtmPrsms` finds the proteoforms containing the tag
      whose mass matches the precursor: the N-terminus a protein N-terminal
      form of `-n` (`allowMod` on the protein; s = the form's `getProtShift`
      plus a sum of 0 to `-S` variable PTM masses within the first peak's
      tolerance) or a truncation start i >= 1 (s = -prefix[i] + such a sum,
      plus a Position N-term PTM mass with one PTM fewer and i >= 2), the
      C-terminus the end e with prefix[e + 1] - prefix[start] + water +
      the form's `getPepShift` or the N-term PTM + a sum of PTM masses = P
      within `compPrecErrorTole`. One simple PrSM per (start, end, form,
      N-term PTM): the range, `n_term_shifts` {the form's shift} (0 for
      truncations), `c_term_shifts` {0}, the truncation N-terminal
      modification shift when used, in the file of its type (complete /
      prefix for forms, suffix / internal for truncations). `VarPtmSlowMatch`
      exits on a form shift that is not an allowed `ProtMod` and reads
      `getCTermShifts()[0]`, so both lists must be set.
    - `multi` and `topmg`: the protein alone (`PtmSlowMatch` and the graph
      alignment find their own diagonals).
- **Files and joining.** The stage's mass filter runs as before; the tag
  filter (`tag_filter_processor::process`, a thread pool over the spectra,
  `OrderedTaskRunner` so the files are in msalign order) writes
  `toppic_tag_filter_<stage>[_<TYPE>]` (`toppic_topmg_tag_filter` for the
  hosted stage); `tag_filter_merge::process` writes `<filter>_tag[_<TYPE>]`,
  the name the search then reads: the filter's candidates first, unchanged,
  then the new tag candidates of the spectrum (`one`: a protein listed with
  the same precursor variant gets the tag shifts added to its lists;
  `var`: the same protein, range, N-terminal shift and truncation
  modification is not repeated; `multi`: the same protein and precursor
  variant is not repeated). For `topmg` the tag file is one more input of
  the `SimplePrsmStrMerge` that builds `toppic_topmg_graph_filter` (the
  best 20 x 5 instead of 20 x 4), so the graph alignment's score-10 rule
  (`SimplePrsmFilter`) applies to the tag candidates too. The multi stage's
  tag filter skips the identified spectra (`skip_spec_ids_`), the hosted
  stage's scores the spectra TopMG searches (`search_spec_ids_`);
  `cleanToppicDir` deletes the files. The console prints "Tag filtering
  (<stage>) - started / finished", the index size, a summary (spectra,
  tags, occurrences, candidates) and the joined counts.
- **Without the option** the run is main d56ebc44's (the first 500 spectra
  of 480_F3_02, setting d: every file; the toppic manual is unchanged, the
  options being hidden).
- **Results on 480_F3_02** (toppic_test experiment 150, report 0150; the
  branch's TopFD, setting d `-s 2 -b var_mods.txt`, 1% FDR, 16 threads; base
  3,374 PrSMs / 2,004 proteoforms / 440 proteins in 19:28): `--tag-filter
  var` 3,376 / 2,006 / 441 (+44 s; 1,578 proteoform candidates for 1,197
  spectra, 1.0% decoys), `one` 3,386 / 2,012 / 441 (+4:09: 24,774 proteins
  for 2,649 spectra, 33% decoys, the one-shift search 120 s against 51 s,
  and a one-shift decoy at E 4.2e-4 lowered the multiple-shift stage's skip
  threshold from 0.0215, so it searched 1,337 spectra instead of 1,112),
  `multi` 3,382 / 2,013 / 441 (+1:02), `var,one,multi` 3,389 / 2,015 / 443
  (+6:11): 23 spectra gained, 17 of them through proteins only the tag
  filter listed (HMGN3 59-99 with 38 fragments at E 1.3e-12, PFD1 42-66,
  TM263 1-37, RL12 124-165, ...), 8 lost at the cutoff (2.41 against 2.50),
  217 PrSMs of the same protein re-explained (112 better E-values, 57
  worse; RL38 2-70 from E 5.2e-5 to 6.2e-46). The hosted TopMG search
  (`--topmg topmg-only`, oxidation, phospho, N-terminal acetyl): 2,996 /
  1,764 / 407 against 2,979 / 1,751 / 406, 17 gained and none lost, but
  35:04 against 16:40, the graph alignment taking 1,482 s against 420 s on
  the 16,537 tag candidates above its score 10 (whole, long proteins). A
  tag filter stage costs 50 s on 4,085 spectra (about 180 tags, 85 k
  occurrences per spectrum). The tag filter lists the protein of about half
  of base's final PrSMs (53% var, 54% one, 48% multi, 60% TopMG; three
  quarters of those with 30 or more matched fragments, a tenth of those
  with 10-19): tags need five consecutive matched cleavage sites, so the
  filter complements the mass filters on strong spectra whose protein the
  index misses and cannot replace them.

## Second alignment of the shift searches (`search/diag/second_align_util`, `--disable-second-align`, issue #614)

Branch `second_align` (toppic_test experiment 144, report 0144 and its
amendment 1; merged into main on 2026-10-08 in the form below). On by
default; toppic's hidden `--disable-second-align` (argument key
`secondAlign` = "false"; `PrsmPara::useSecondAlign`; the parameter block
prints "Second alignment at the refined precursor: false" only then)
keeps the first alignment's PrSM, main's behaviour before the merge.
Facts to keep straight:

- **Why.** The one-shift and multiple-shift alignments (PsAlign,
  PsAlignMultiType) match the precursor-dependent ions with the fragment
  tolerance plus the precursor's (the previous section), and the PrSM's
  annotation matches them with the fragment tolerance at the refined
  precursor, so the dynamic program counted ions the annotation did not
  and could choose a diagonal set on their strength (report 0135 section
  7: "A second alignment on the refined precursor with the strict
  tolerance would make the alignment and the annotation agree").
- **What it does** (`PsAlign::geneResult`, `PsAlignMultiType::geneResult`).
  The first pass is the former `geneResult`, factored into
  `second_align_util::genePrsm` (the precursor refined by the matched
  ions to P', the extended spectra at P', the ranges refined, the PrSM).
  Then, for every call with a first alignment: the PRM peaks are rebuilt
  at P' (`genePrmPeaks`: `prm_ms_factory::geneMsSixPtrVec` with the
  refined precursor, so the REVERSED peaks are P' - (m - C shift), plus
  the zero and precursor peaks sized on P'; every fragment peak's
  N-strict-C-relax and N-relax-C-strict tolerances set to its strict
  tolerance ppo x m); the first alignment's headers are moved to P'
  (`moveHeaders`: an N-anchored header keeps its N-terminal shift, a
  C-anchored one keeps its C-terminal shift, so its N-terminal shift
  moves by P' - P; C-terminal shifts and truncation positions from
  `DiagHeader::initHeader` at P'); the diagonals are rebuilt
  (`geneDiagonals`, headers without a pair dropped); the same dynamic
  program runs on them (PsAlign for the type of the last `compute()`,
  which it now remembers; PsAlignMultiType for the types of the last
  `compute()` together, kept by P' in `second_align_ptrs_` so the types
  whose first results refine to the same mass share one second
  alignment); the second result for the same type and shift number goes
  through `genePrsm` with `keep_prec` (the refinement's delta range is
  [0, 0], so P' is kept and the headers keep their shifts), the range
  refinement and the PrSM. That PrSM is returned; with no second
  alignment of positive score the first PrSM is kept, so no candidate is
  lost. The candidates are ranked by the PrSM's matched fragments as
  before, and the E-values are computed from the PrSM afterwards (the
  generating function is unchanged). The zero-shift search, the variable
  PTM search and TopMG's graph alignment are untouched. The precursor
  never leaves the window of the first refinement (+-T around the
  observed precursor).
- **`refinePrecursorByMatchedIons` takes an optional delta range**
  (`delta_min` / `delta_max`, infinite by default), which cuts its
  [-T, T]; only the second pass uses it, with [0, 0].
- **Research log.** `TOPPIC_SECOND_ALIGN_LOG=<file>` appends one
  tab-separated line per `geneResult` call with a first alignment: stage
  (`one` / `multi`), spectrum id, type, shift number, protein, P, the
  first score, P', the first PrSM (start, end, matched theoretical
  fragments, shifts as left:right:mass), the second score (`-` without a
  second alignment), P' again, the second PrSM, and `second` / `first`
  for the PrSM kept. The alignment's score is the sum of the diagonal
  pairs' scores (one per spectrum and ion type per matched cleavage) with
  the zero and precursor peaks counted with the group's spectrum number
  each, so score - 2 against the PrSM's matched fragments measures the
  agreement (for the multiple-shift search add the large-shift penalty
  back, 3 per shift above 300 Da).
- **Variants tested and dropped** (experiment 144, branch commits
  3dc127c4-26f7dbe3, kept in the branch's history): a second refinement
  by the matched ions after the second alignment (+-T around P', so up to
  2T from the observed precursor: 14% of the one-shift candidates and 68
  reported shift PrSMs of 480_F3_02 ended beyond 10 ppm, and one such
  move lost two ATIF1 oxidation PrSMs through the cluster protein
  correction), the same clamped to the observed precursor's window, and
  a second alignment with the first alignment's relaxed tolerances
  (which left the alignment and the annotation apart: the tolerance is
  the gap, not the precursor).
- **Results on 480_F3_02** (setting d, experiment 131's TopFD output,
  against main 7b5112c2): the second alignment's score is within 1 of the
  PrSM's matched fragments for 98.4% of the 265,195 one-shift and 98.3%
  of the 41,712 multiple-shift candidate alignments (the first alignment:
  81.9% and 69.3%). It changes the positions or shifts of 14% and 19% of
  the candidate PrSMs (the shift positions and diagonal sets; 695
  one-shift candidates gain matched fragments, mostly +1, and 297 lose).
  1% FDR: 3,356 PrSMs / 1,996 proteoforms / 439 proteins against 3,355 /
  1,997 / 438, the same E-value cutoff (2.49) and decoys: within the data
  set's noise, because the candidates that reach 1% FDR are chosen by the
  annotation's fragments and the E-value either way, and the relaxed
  window inflated targets and decoys alike. The one-shift search takes
  38 s instead of 23 s; the multiple-shift search's 256 s are unchanged.
  The pushed main 249433b1 (this, the penalty fix below, `filter_pairs`,
  `topmg_prec_matched_ions` and `local_align`) gives 3,375 / 2,004 / 440
  against the pre-merge main c7f13d66's 3,372 / 2,008 / 439, the same
  cutoff (report 0144 amendment 2).

## Large-shift penalty of the shift alignments (`search/oneptmsearch/ps_align`, `search/ptmsearch/ps_align_multi_type`)

Branch `ps_align_penalty` (toppic_test experiment 145, report 0145; merged
into main on 2026-10-08 as 0bca905d, after `second_align` e61e339b; on the
development subset the merged main with `--disable-second-align` differs
from main c7f13d66 only in the one-shift stage files and what follows them,
as the branch did from its base). `PsAlignPara` fixes `align_large_shift_thresh_` = 300 Da and
`align_large_shift_panelty_` = 3: a shift between two diagonals above the
threshold costs 3 points in the dynamic program. Facts to keep straight:

- **The two programs applied it differently until this branch.**
  `PsAlignMultiType::dp` (the multiple-shift search since 4f845a6f)
  compares the shift predecessor's score minus the penalty with the
  diagonal and truncation alternatives and stores it, so the penalty
  propagates and the end score is 3 lower per large shift. `PsAlign`
  (the one-shift search; the multiple-shift search with
  `--disable-multi-type-align`; the original TopPIC code) used the
  penalty in `getShiftPre` only to choose which diagonal to shift from,
  and `dp` then took the chosen predecessor's score without it: a large
  shift was never discouraged against continuing without a shift, and
  the stored and final scores carried no penalty. `ps_align_multi_type.hpp`
  describes the two as the same recursion; they were not (toppic_test
  report 0144 section 7, found by comparing the alignment scores with
  the PrSMs' matched fragments).
- **This branch makes PsAlign match PsAlignMultiType**: `getShiftPre`
  returns the chosen predecessor's penalized score, and `dp` compares and
  stores it. The predecessor choice is unchanged; what changes is that a
  path through a large shift now competes with the no-shift and
  truncation alternatives 3 points down, and the end score reflects it.
- **Results on 480_F3_02** (toppic_test experiment 145, report 0145;
  setting d, experiment 131's TopFD output, against main 7b5112c2): 3,354
  PrSMs / 1,997 proteoforms / 437 proteins against 3,355 / 1,997 / 438,
  the same E-value cutoff. 825 of the 14,670 best one-shift candidate
  alignments change, 672 of them a shift above 300 Da replaced by a
  smaller one (the same matched fragments in 593, 1-3 fewer in 206, more
  in 26; targets and decoys alike); the multiple-shift stage files are
  unchanged. 31 of 3,353 reported PrSMs get another proteoform, mostly a
  large shift rewritten as a smaller one with the same or better
  fragments and E-value; 2 MED9 PrSMs are lost when the feature's YBOX1
  +14 Da PrSM takes the cluster (the penalty took 2905's -421 Da
  alignment with 19 fragments down to a -49 Da one with 18) and SEM7A at
  the protein cutoff; 1 gained. No time cost.

## Tolerance of the complementary PRM peaks in the E-values (`stat/tdgf/prob_peak`, `--evalue-complement-tolerance`, issue #616)

Branch `evalue_strict_complement` (toppic_test experiment 146, report 0146;
from main 90d744d7; merged into main on 2026-10-08 as a hidden option with
the relaxed tolerance, main's behaviour, as the default). Facts to keep
straight:

- **The three tolerances of a PRM peak** (`prm_ms_factory::addSixMasses`,
  the fragment tolerance ppm): an ORIGINAL peak (a b-ion mass m) has
  strict = N-strict-C-relax = ppm x m and N-relax-C-strict = ppm x (m + P);
  a REVERSED peak (the complementary mass P - m of a y-ion, P the
  precursor the PRM spectrum was built at) has strict = N-relax-C-strict
  = ppm x m and N-strict-C-relax = ppm x (m + P), the precursor's error
  added because P - m moves with P. The searches use the relaxed
  tolerances for the precursor-dependent ions (`compDiagPair`), the
  annotation and the second alignment the strict one at the refined
  precursor.
- **The generating function before the branch** (`ProbPeak`, the "strict"
  flag false in toppic; the separate mode with the PRM spectra at the
  PrSM's precursor, `CompPValueArray::compSingleExpectedValue`, is not
  used) took every peak's N-strict-C-relax tolerance: the strict window
  for an ORIGINAL peak and the relaxed one, about twice as wide, for a
  REVERSED peak, on the six-mass PRM spectra of the spectrum set at the
  observed precursor P (or P +- 1 isotope for the PrSMs found there).
- **The option**: toppic's hidden `--evalue-complement-tolerance
  <relax|strict|separate>` (argument `evalueComplementTolerance`, default
  `relax`; anything else is rejected). With `strict`, `ProbPeak` takes
  the strict tolerance for a REVERSED peak too
  (`PrsmPara::useEvalueStrictComplement`, passed to `CompProbValue` by
  `CompPValueArray`; the cluster protein correction's recomputation goes
  through the same `PrsmPara`); the zero and precursor peaks are ORIGINAL
  and unchanged, so the precursor peak's width, which normalizes the
  shift E-values (`compPrecProbs`), stays the precursor tolerance. The
  parameter block prints "E-value tolerance of complementary ions:
  strict" (or separate) only when set. topmg has no option (its own
  E-values are the MCMC's; `--topmg` in toppic follows the setting).
  Without the option the run is main's before the merge byte for byte
  (development subset: 63 files but the time stamps).
  `separate` (6f612533) selects the code's separate mode instead
  (`PrsmPara::useEvalueSeparate`, `EValueProcessor::process(true)` in
  both E-value stages and the cluster protein correction): one
  generating function per PrSM on PRM spectra rebuilt at the PrSM's
  adjusted precursor, every peak strict, the proteoform mass as the
  precursor. It costs nothing (the E-value stage's time is the table
  setup per spectrum, 176 s on 480 either way).
- **The generating function's resolution** (report 0146 section 7):
  `TdgfMng::convert_ratio_` is 27.4335 units per Da (0.0365 Da; the
  searches use 274.335), residue masses are rounded to integers (proline
  +0.50 unit), and a window is ceil(tolerance x 27.4335) units, so a
  strict window of 10 ppm x m is one unit for every fragment below 3.65
  kDa and a relaxed one, 10 ppm x (m + P), one unit below 3.65 kDa and
  up to six at 20 kDa. The strict setting therefore changes nothing for
  precursors below 3.6 kDa and most for the heaviest, and a PrSM whose
  ions form a long ladder on a proline- or D/E-rich stretch can get a
  probability of 0 ("Invalid e value", the placeholder 1e300) because no
  integer path hits all its one-unit windows (NPM scans 3786 and 3617,
  FNBP4 6285 on 480; the separate mode has them too).
- **Results on 480_F3_02** (setting d, experiment 131's TopFD output,
  against main 8688fd08): strict 3,374 PrSMs / 1,994 proteoforms / 441
  proteins, separate 3,373 / 1,997 / 444, base 3,375 / 2,004 / 440.
  Every E-value falls (targets median x 10^-0.45, decoys x 10^-0.31,
  more for heavier precursors and higher scores) and the cutoffs with
  them (PrSM 0.84 against 2.50); 8 spectra report a one-shift form of
  the same protein instead of main's plain form, the multiple-shift
  stage skips 41 more spectra because its skip threshold (the smallest
  first-pass decoy E-value, experiment 131) falls from 0.021 to 0.009
  and 7 two-shift PrSMs are lost that way (ROA2 5019 with 39 fragments,
  JUPI1 4202 with 40). Not merged; the recommendation is `relax` as the
  default if the option is merged, and a finer generating-function
  resolution before the strict tolerance is tried again.

## Range of a variable PTM shift (`seq/mass_shift`, `search/varptmsearch/var_ptm_align`)

Branch `var_ptm_narrow_range` (toppic_test experiment 122, report 0122;
merged into main on 2026-10-05). Facts to keep straight:

- **The ranges.** `VarPtmAlign::geneShiftVec` takes the variable PTM shifts
  from the gaps between the refined diagonals
  (`diag_header_util::getDiagonalMassChanges`) and narrows each alteration
  to the residues of its gap that can carry its PTM. It now also calls
  `MassShift::narrowRangeToAlters`, which sets the shift's range to the
  union of its alterations' ranges (alterations covering no residue are
  skipped; with none left the range stays as it is). A shift with several
  alterations (a diagonal between two PTMs was empty) covers all of them.
- **All residues of the PTM** (amendment 1). An alteration covers every
  residue that one of the modifications of its PTM modifies:
  `VarPtmSearchMng::ptm_mod_ptr_vec_2d_` keeps them in one list per PTM
  (grouped by PTM, not by mass), `getPtmModPtrVec(mod)` returns it, and
  `geneShiftVec` tests each. `Deamidated,0.984016,NQ,...` gives Deamidated
  on N and on Q; the alignment puts a tied PTM on the first N or Q of its
  stretch (the `>` tie rule of `VarPtmAlign::dp`) and the backtrack takes
  that residue's modification, so testing only it (as before) left the
  other type out of the range (FOXK1 4009 `AVT(Q)[Deamidated]NGK`, now
  `AVT(QN)[Deamidated]GK`). On 480_F3_02 (ox + deamidation) 3,247 first-pass
  ranges and 8 reported ones change; search fragments, E-values and
  identifications do not, and oxidation-only runs are identical.
- **What the range decides.** `Proteoform::getSeqSegmentPtrVec` generates
  no ion inside a shift's range, and the proteoform string and the
  "Variable PTMs" column show it (`P(M)[Oxidation]`, `Oxidation:[44]`
  instead of `(PM)[Oxidation]`, `Oxidation:[43-44]`).
- **No ion of the search spectrum is gained.** The range refinement trims
  each diagonal to its matched ions in the same extended spectrum, so the
  ions the narrowing adds were already unmatched there. On 480_F3_02 the
  first-pass PrSMs keep their counts and E-values (10,791 / 11,854 variable
  PTM PrSMs, a Python copy of `findPairs` predicts every count).
- **Fragments are gained after the search.** The post mass matching looks
  for every unmatched theoretical fragment of the best PrSM in TopFD's raw
  peaks, and the recount uses every mass. With more theoretical ions both
  find more: of the best PrSMs of every spectrum, 203 / 231 gain 482 / 541
  final fragments (88-89% of them masses the post mass matching added);
  none loses. E-values come from the search, so they do not change, except
  for a PrSM that the cluster protein correction gives another spectrum's
  proteoform (the correction recomputes its E-value): 1 of 4,075.
- **Identifications are unchanged** (oxidation 3,336 / 1,998 / 437,
  oxidation + deamidation 3,336 / 1,992 / 436), as are the decoys and the
  run time. Of the reported variable PTM shifts, 438 / 480 (ox) and 543 /
  592 (ox + deamidation) are on one residue (main: 143 and 208); every
  other one covers two or more residues that can carry its PTMs.
- **Determined sites** (branch `var_ptm_split_sites`, toppic_test
  experiment 124, report 0124; merged into main on 2026-10-05). A shift holds several PTMs when the
  diagonal between them had no matched fragment. When each of its PTMs has
  as many candidate residues as copies (and no residue is a candidate of
  two PTMs), the sites are known: `VarPtmAlign::splitDeterminedShift` makes
  one shift per PTM on its residue, with that residue's modification
  (`(MAQM)[Oxidation;Oxidation]` -> `(M)[Oxidation]AQ(M)[Oxidation]`,
  `Oxidation:[33];Oxidation:[36]`). Any other shift is left as it is. The
  fragments between the sites cannot match in the search spectrum (that
  diagonal was empty), so on 480_F3_02 no first-pass count or E-value
  changes (1,589 / 2,451 first-pass PrSMs split); 11 reported PrSMs per
  run change their notation, 3-4 of them gain a final fragment.

## Deamidation or the -1 isotope offset (`prsm/prsm_deamid_selector`)

Branch `deamid_prec_error` (toppic_test experiment 123, report 0123;
merged into main on 2026-10-05). **Off by default**; the hidden toppic option `--deamid-selector`
(key `deamidSelector`, default `"false"`; printed in the parameter block
only when set) turns it on. Then, with a deamidation in the variable PTM
file (a PTM with Unimod id 7 or a monoisotopic mass within 0.001 Da of
0.984016 Da; `getDeamidPtmNames`), `TopPIC_identify` runs
`prsm_deamid_selector::process` between the E-value computation and the top
PrSM selection: it reads `toppic_evalue` and writes `toppic_evalue_deamid`,
which `prsm_top_selector` then reads. Facts to keep straight:

- **Why.** The unmodified proteoform matched at the -1 offset (a zero-shift
  PrSM whose adjusted precursor mass is P - 1.00228) and the deamidated one
  at P (a variable PTM PrSM) explain a spectrum with masses 18.3 mDa apart.
  Since experiment 117 their E-values come from generating functions on the
  PRM spectra of different precursor masses, and the offset match often
  wins with fewer fragments.
- **The pairs.** Per spectrum (all its precursors, as the top PrSM selection
  groups them) and per protein, start and end: the zero-PTM PrSMs
  (`unexpected_ptm_num` and `variable_ptm_num` 0) and the deamidation PrSMs
  (no unexpected shift, a `Variable` alteration whose PTM abbreviation is a
  deamidation, read with `MassShift(XmlDOMElement)`). The best of each kind
  by E-value is compared. The pairs are zero-PTM PrSMs at the -1 offset in
  practice: no other offset gives the two masses within the tolerance.
- **The rule:** matched fragments (`match_theo_frag_num`, the E-value's
  score) 5 or more apart decide. Otherwise the smaller |precursor error|
  decides: the label-free precursor mass + k x 1.00228 - the adjusted
  precursor mass (the proteoform mass in both searches), k = round((adjusted
  - precursor) / 1.00228). A tie keeps the better E-value. Every PrSM of the
  losing kind in the group is removed, so that another deamidation site
  cannot take the loser's place. E-values are not changed.
- **The loser may have been the best PrSM, the winner need not be.** When
  the zero-PTM PrSM wins with a worse E-value than the deamidation PrSM, the
  top PrSM selection can pick a third PrSM (on 480: one-shift PrSMs of the
  same protein).
- **Results on 480_F3_02** (`-s 2`, oxidation and deamidation, against main
  7c8e00e6; TopFD 47b9fc57 / ead67ddc outputs; run with the step on, before
  the option made it opt-in):
    - The step: 539 / 544 pairs, half of them decoys, 620 / 621 PrSMs
      removed, 3.4 s. A Python copy of the rule predicts
      `toppic_evalue_deamid` exactly.
    - 3,335 / 1,994 / 436 and 3,338 / 1,997 / 436 at 1% FDR (main 3,335 /
      1,993 / 436 and 3,338 / 1,996 / 436), the same cutoff and decoys.
    - It recovers report 0117's five HMGN1/HMGN2 deamidations (2488, 2752,
      2833, 2838, 2974) on both inputs and keeps 4679's offset match on D.
    - It loses 7 other deamidations per input: PTMS 33-80 (5,553 Da; five
      scans on S, six on D) and PTMA (12 kDa) by precursor error, LMNA 3734
      (S) by fragments; one-shift PrSMs are reported instead. Deamidated
      PrSMs at 1% FDR 94 -> 92 and 102 -> 100.
- **The cause of the PTMS losses:** TopFD's precursor masses are
  systematically high by about +1.3 to +1.6 ppm (+7 mDa at 5.5 kDa, +17.5
  mDa above 10 kDa; the median of the confident zero-shift PrSMs), so the
  deamidated explanation (+10 to +13 mDa) looks worse than the -1 offset
  (-5 to -8 mDa). Measured from the run's median ppm error instead of zero,
  the deamidation is closer in all five (report 0123's simulation; not
  implemented).

## Truncation N-terminal modifications (`-b`, Position `N-term`)

A variable PTM file line with Position `N-term` is a variable PTM of the
proteoform N-terminus (the alpha-amine of its first residue), kept apart
from the side-chain PTMs of that residue. Facts to keep straight:

- **Representation.** `Proteoform::trunc_n_term_mod_ptr_` (default
  `ModBase::getNTermNoneModPtr()`) records it, and its mass is the residue
  sequence's `n_mod_ptr_` (sequence mass, `BpSpec` prefix masses), as for an
  acetylated protein N-terminal form; it is no mass shift.
  `setTruncNTermModPtr` accepts it only on a proteoform with the NONE
  protein N-terminal form starting at protein residue 3 or later (0-based
  start >= 2; residues 1-2 are the protein N-terminal forms') whose first
  amino acid is the mod's. `proteoform_factory::geneTruncNTermModProteoform`
  makes the variant, `geneSubProteoform` keeps it when the first residue
  stays, the PrSM XML has `<trunc_n_term_mod>` only when set.
- **Output.** `[Name]-` prefix of the proteoform string, outside any
  annotation of residue 1; `Name:[N-term]` first in `Variable PTMs`
  (`getVarPtmStr`); `getVarPtmNum` (`# Variable PTMs`, XML
  `variable_ptm_num`, FDR groups) counts it; one `prsm_mass_shift` row of
  type `N-terminal variable`.
- **Filter.** One more round per N-term mass m (spectrum built as for an
  N-terminal label of m), only on diagonal-index rows of a fitting start
  (suffix and internal types; the N-terminal form index is not scored), up
  to `-S` - 1 other PTMs, each round its own `SimplePrsmMemMerge` group. A
  candidate carries the masses to try in `SimplePrsm::trunc_n_term_mod_
  shifts_` (`{0}` = none, not written to the file then); a candidate of
  several rounds is listed once with the union.
- **Search.** The alignment knows no N-term shift: the processor searches
  the candidate for each listed mass, 0 as it is, otherwise the variants of
  the N-term mods of that mass, with at most `-S` - 1 side-chain PTMs
  (`VarPtmAlign::max_ptm_num_`). So residue 1 can carry an N-terminal and a
  side-chain PTM. When a variant and the candidate as it is score the same
  (one acetylation of K1, on the N-terminus or the side chain), the
  processor's sort (`cmpScoreDecProtIncTruncNTermFirst`) reports the
  N-terminal one; the user chose this over comparing E-values.
- **E-value.** (T + 1)^k x (T_N + 1) for a PrSM with k side-chain PTMs and
  a truncation N-terminal modification, T_N = the number of N-term masses
  (`TdgfMng::n_term_var_ptm_type_num_`, argument `nTermVarPtmTypeNum`); the
  candidate count uses the sequence mass without the N-terminal mod.

## Proteoform classification (`seq/proteoform`, the PrSM tables)

Branch `proteoform_classification` (issue #611; toppic_test experiment 129
and its amendment 1, report 0129; merged into main on 2026-10-06). `Proteoform::getClassification(
prot_hit_num, var_mod_ptrs)` gives the level of the proteoform in the
five-level classification of proteoform identifications (Smith et al., Nat.
Methods 16:939-940, 2019), by the rule the user gave:

|  | one protein | several proteins |
|---|---|---|
| no unexpected shift, variable PTM sites decided | `1` | `2D` |
| no unexpected shift, a variable PTM site undecided | `2A` | `3` |
| unexpected shift(s) | `4` | `5` |

Facts to keep straight:

- **The properties.** An unexpected shift is an alteration of type
  `UNEXPECTED` (`getAlterNum`); protein N-terminal forms, fixed PTMs, the
  truncation N-terminal modification and TopMG's FASTA PTMs (`INPUT`) do
  not count. The protein hits are `SearchFastaMatch::process`'s matches,
  the `# Protein hits` column; 0 counts as one protein. A variable PTM site
  is undecided (`hasUndecidedVarPtmSite`) when a variable mass shift covers
  more than one residue and one of its PTMs does not have exactly as many
  candidate residues as alterations, or two PTMs share one: the candidates
  of a PTM are the residues of the range that one of its modifications in
  `var_mod_ptrs` modifies (Deamidated on N and Q, Phospho on S, T and Y),
  the test of `VarPtmAlign::splitDeterminedShift`. toppic splits a shift
  whose sites are decided (one PTM per residue), so in toppic every
  multi-residue variable shift is undecided; TopMG's shifts span the gaps
  between matched fragments, and most of its ranges hold one candidate
  residue (decided). The proteoform knows neither the database nor the PTM
  file, so both come from the caller.
- **One place computes it.** `prsm_cutoff_selector::process(...,
  var_mod_ptrs)` sets it for every selected PrSM without one, from
  `SearchFastaMatch` on the database file the tables search (built only
  when a PrSM needs it); `TopPIC_report` passes the "anywhere"
  modifications (`mod_util::readAnywhereModTxt`) of `-b` and `-B`,
  `TopMG_post` those of `-i`. The proteoform stores it (`setClassification`,
  `getClassification()` without arguments; "" until set). The proteoform
  and protein cutoffs and the top PrSM filters read the PrSM cutoff's PrSMs
  and keep it, and the tables and the SQLite writer, which read the cutoff
  files, write the stored level (`-` if none).
- **Where it is written.** The last column, `Proteoform classification`, of
  every table of `PrsmMatchTableWriter` (PrSM, proteoform and protein,
  single and with the other matches; toppic, its combined `-c` tables and
  topmg), `proteoform_classification`, the last column of the SQLite tables
  `prsm`, `proteoform` and `protein`, and `<classification>` in the
  `<proteoform>` element of the XML outputs (`_prsm.xml`, `_proteoform.xml`,
  `_protein.xml`), after `<proteo_db_seq>` and the truncation N-terminal
  modification. The column is last so that no existing column moves; the
  rows of a PrSM's other protein matches have no field for it (they end
  before the q-values, 39 of the 45 fields before it). The element is
  written and read only when set, so the PrSM files before the report are
  unchanged.
- **What the rule leaves out.** Levels 2B and 2C never occur. An unexpected
  shift is level 4 whatever its mass, also when it is a known PTM on one
  residue (HMGN1 `(S)[+79.9669]`) or of isotope size. With `-B`, a shift
  that the localization explains by a PTM becomes a `VARIABLE` alteration
  (`stat/local/local_proteoform`), so its proteoform is level 1, 2A, 2D or
  3 (not run). Several proteins can be one gene (LAP2A / LAP2B, both TMPO),
  and that is still 2D. The level is that of the PrSM's proteoform, so the
  PrSMs of one proteoform id can differ (one with an unexpected shift, one
  without), and the proteoform table shows the level of its best PrSM.
- **Results** (toppic_test experiment 129, 480_F3_02, `-f C57 -d -t FDR
  -T FDR`, oxidation): toppic (`-y FDR -s 2`) writes main 4ecb19ea's files
  but the column, the element and the SQLite column (1de20d8b, levels 1 /
  2D / 4 / 5, against main; a2a84e16 against 1de20d8b: the report stage's
  PrSM files equal without `<classification>`, the tables and SQLite
  tables equal but 29 PrSMs, 16 proteoforms and 1 protein moved from 1 to
  2A), in the same time (16:48 against 16:42; the PrSM cutoff 1.4 s
  against 0.3 s, the protein hit index). The level agrees with a Python
  copy of the rule in every PrSM row of toppic's and topmg's tables, and
  the XML element and the SQLite column with the tables. Proteoforms of
  levels 1 / 2A / 2D / 3 / 4 / 5 at 1% FDR: toppic 1,050 / 16 / 70 / 0 /
  812 / 50; topmg (one unexpected shift by default) 1,007 / 15 / 68 / 1 /
  21 / 1, where 210 of its 228 range shifts hold as many candidate
  residues as PTMs (decided).

## PTM localization by a variable PTM search (`stat/local_align`, the default; `--local-original` for `stat/local`)

Branch `local_align` (toppic_test experiment 140, report 0140; from main
6eb40e97, merged into main as 927e2da3 on 2026-10-08). Since 2026-10-08
`LocalAlignProcessor` is **the default** `-B` localization (argument
`localAlign`, "true"); the hidden `--local-original` selects the MIScore
candidate search of `stat/local` (`LocalProcessor`, the method before),
and the hidden `--local-align` (the switch that selected the new method
while it was optional) is still accepted. The parameter block prints
"Localization method: Variable PTM search (local_align)" or "MIScore
candidate search (local, --local-original)". Without `-B` neither runs.
What the method does, and how it differs from `stat/local`:

- **Explanations of a shift.** For an unexpected shift of mass m and the
  precursor tolerance e (`compStrictErrorTole(adjusted precursor mass)`),
  the explanations are the single PTMs X of the `-B` file with
  |mass(X) - m| < e and a modifiable site, and the unordered pairs {X, Y}
  (X = Y allowed) with |mass(X) + mass(Y) - m| < e and two sites in one of
  the two orders. Unlike `local_util::getPtmPtrVecByMass`, a match one
  isotope mass away does **not** count by default (switch `iso`, 0 since
  experiment 140's amendment 1: with it on, a +12.98 Da shift of a PrSM
  with a +1.00 Da second shift became a Methyl at -1 Da, i.e. a 1 Da
  precursor error turned into a PTM). The shift is replaced only when
  there is **exactly one** explanation; with none or several it stays
  (logged "none" / "ambiguous"). Several explanations are first reduced
  by the **tie-break rules** of `local_align_util::breakTies` (switch
  `tie`, 1; experiment 140's amendment 5), applied in this order while
  more than one is left: (1) fewer PTMs (Dimethyl over Methyl+Methyl,
  Trimethyl over Methyl+Dimethyl, a PTM on the N-terminal acetylation
  form over Acetyl + the PTM on the plain form), (2) the N-terminal
  acetylation form over an Acetyl of `-B` on the plain form of the same
  start (terminal candidates); (3) **Acetyl vs Trimethyl** (0.036 Da
  apart, within 10 ppm above 3.6 kDa): when the explanations left are on
  the same proteoform and the same up to Acetyl <-> Trimethyl, both are
  localized and the better match score wins (then the higher confidence,
  then Acetyl; `LocalAlignProcessor::choose`, logged
  `site_scores(Acetyl:31/Trimethyl:33)`) - a first version of amendment
  5 preferred Acetyl outright, which gave calmodulin's K115
  trimethylation (758, 88.5%) an Acetyl at 48.9%. What is left after the
  rules must be one explanation. Without the rules (`tie=0`) Dimethyl is
  never chosen when Methyl is listed too, nor Acetyl when Trimethyl is.
  The PTMs are recognized by the abbreviations `Acetyl` and `Trimethyl`
  of the `-B` file. The log's `rule` column names the rules that decided
  and `all_explanations` the set before them.
- **No filters.** `LocalProcessor` keeps a localization only when the
  matched fragment number stays above the original's minus 5 and above
  0.67 x the original (`DESC_MATCH_LIMIT_`, `desc_ratio_`), when a site
  passes the MIScore threshold `-H`, and (two shifts) when the localized
  PTMs overlap the original shifts; `LocalAlignProcessor` applies none of
  these. When no site reaches `-H`, the PTM is reported on the range from
  the first to the last site with a positive score, with the **best
  site's score** as the confidence (`local_align_util::filterScores`;
  `LocalAnno::getConf`, the `<confidence>` of the PrSM XML). Above the
  threshold the confidence is `scrFilter`'s, the sum of the scores of the
  reported range; the sum over the fallback range would be 1.0 whatever
  the spread, so it was replaced by the best site's score in experiment
  140's amendment 1 (a Methyl over 18 K/R of H2B at 5.5% each: 0.055, not
  1.0). The per-site scores are in the MIScore column either way.
- **The proteoform's termini.** `stat/local` scores the candidate forms
  of `local_proteoform::getAllCandidateForms` (N- and C-terminal
  truncations by up to 5 residues, the N-terminal acetylation forms) and
  takes the best; `local_align` first explains the shift value of the
  proteoform as it is. When that gives no explanation and the shift has
  at most `term_max_matched` (3; 2 until amendment 4) matched cleavage
  sites on its N-terminal (C-terminal) side (`compMatchedSites` of the
  original proteoform), the candidate proteoforms whose start (end) is
  moved by up to `term_range` (5) residues are tried
  (`LocalAlignProcessor::processTerminal`, the forms of
  `local_proteoform::getCandidateForm`, so a new start of 0 or 1 brings
  the N-terminal acetylation variants; the piece keeps at least one
  residue, and a truncation may cut into the shift's range - until
  amendment 4 the shift's residues had to stay inside, which skipped the
  truncations that explain a shift sitting on the terminal residues, 31
  of main's moved-start localizations on 480). The other N-terminal forms
  of the **unmoved** start (the acetylation dropped, or added: main's base
  candidates include them, and +37.95 Da + 42.01 is a Phospho) belong to
  the **standard attempt** (`processPiece`, `n_forms`): when the piece
  holds the N-terminus and the shift has at most `term_max_matched`
  matched sites on its N-terminal side, their explanations join the pool
  of the proteoform's own before the tie-break rules, so that a PTM on the
  acetylation form competes with Acetyl + the PTM on the plain form
  (amendments 4 and 5; until 8476656b they were terminal candidates, tried
  only after an "ambiguous" standard attempt had already stopped). Each
  candidate's unexplained mass is the adjusted
  precursor mass minus the candidate's mass, and the one-explanation rule
  counts over all (candidate, explanation) pairs ("term_localized",
  "term_none", "term_ambiguous" in the log, the explanations tagged
  `@start<n>` / `@end<n>` and the N-terminal form). The result proteoform
  takes the moved terminus (the base form from `getCandidateForm` with the
  kept shifts, the new shifts moved with `shiftPositions`); a two-shift
  PrSM can move its start for the first piece and its end for the second.
  Switch `term` (1; 0 = the proteoform as it is only). Added in
  experiment 140's amendment 3.
- **The site tables** (`compOnePtmSTable`, `compTwoPtmMTable` in
  `local_align_util.cpp`) are `local_util`'s with the known shifts of the
  proteoform (the variable PTMs, and the other unexpected shift of a
  piece) in the theoretical prefix and suffix masses, a cleavage site
  inside a known shift's range giving no fragment (as
  `Proteoform::getSeqSegmentPtrVec` does for the matching);
  `local_util::compOnePtmSTable` ignores the known shifts. The protein
  N-terminal modification of a proteoform read from the PrSM file (the
  acetylation) is the residue sequence's N-terminal mod
  (`ResidueSeq::getNModPtr`), not a residue's mass nor a shift that
  `Proteoform::getMass` adds; `makePiece` adds it to the prefix masses of
  a piece starting at the N-terminus (as `BpSpec` does) and to the N
  offset of any other piece - missed until experiment 140's amendment 7,
  so the N-terminal fragments of acetylated proteoforms never matched
  (13 of 39 acetylated PrSMs localized by both main and the branch had
  different site scores, none of the 194 plain ones). The fragment
  match is a binary search of the sorted extended masses per theoretical
  mass (not `compMTable`'s two-pointer walk, which needs monotonic
  theoretical masses, broken by a negative known shift). The scores are
  `prob_ratio_ ^ (score - best score)` before the normalization, so that
  a score above ~230 does not overflow `pow`.
- **Two PTMs for one shift**: both orders are scored with the DP of
  `LocalProcessor::compTwoPtmScr` and the better order is localized by
  the break-point method of `twoPtmLocalize` (ported, with the site pairs'
  probabilities summed directly). The two VARIABLE shifts take the ranges
  of their sites above the threshold, `[bgn_1, end_1]` and
  `[bgn_2, end_2]`. `local_proteoform::createProteoformPtr` builds its
  alterations with the break point as one end, but `Alter::setLocalAnno`
  resets an alteration's range to its annotation's `[bgn, end + 1)`, so
  main's ranges are the same tight ranges and match the score vectors
  that `Proteoform::getMIScore` and the ProForma writer index from the
  alteration's left position. Experiment 140's amendment 2 first read
  this as a latent out-of-range read and padded the vectors to the
  break-point ranges (4cdb0069), which emptied the MIScore of the second
  PTM of every two-PTM PrSM (13 on 480); the padding was reverted and the
  comment in `createProteoformPtr` states the rule.
- **Two unexpected shifts.** `stat/local` explains the sum of the two
  shifts by one ordered pair over the whole proteoform; `local_align`
  splits the proteoform at a cleavage site b between the shifts (first
  shift's right position <= b <= second's left position): the middle one
  (by index) of the sites in that range matched by a fragment of the
  spectra (`compMatchedSites`, the original proteoform with all its
  shifts), or the middle of the range when none is matched. The piece
  before b holds the first shift, the piece from b the second
  (`makePiece`; a known shift belongs to the piece of its left position
  and is cut at the piece's end), and each piece is explained on its own
  with the residue masses and shifts of the other piece as its fragment
  offsets (`Piece::n_offset`, `c_offset`), so that the pieces are matched
  against the real spectra. A PrSM can come out with one shift explained
  and the other kept (switch `split`, 0 = two-shift PrSMs left alone).
- **The result PrSM** keeps the proteoform's residues, N-terminal form,
  truncation N-terminal modification, start / end and cluster ids; the
  explained shifts are replaced by the VARIABLE shifts, the PrSM takes the
  proteoform with `setProteoformPtr` (matched numbers recounted) and the
  adjusted precursor mass becomes the proteoform's mass, as `stat/local`
  does. A one-shift PrSM whose shift is within [-1, 1] Da + e is left
  alone (`stat/local`'s rule).
- **Alternative explanations** (`Proteoform::alternative_mass_shift_list_`,
  `MassShiftPtrVec2D`; toppic_test experiment 140, amendment 8). When a
  shift has several explanations, the other candidates on the same
  proteoform as the chosen one (the proteoform as it is when none is
  chosen: the shift stays, the alternatives are all its explanations) are
  localized too (`LocalAlignProcessor::choose`, `PieceResult::alternatives`)
  and the result proteoform gets one complete mass shift list per
  alternative (`setAlternatives`): its shifts other than the unexpected
  ones and the variable PTMs the localization added, the alternative's
  variable PTMs, and any other unexpected shift as it is - so for a
  two-shift PrSM each list explains one shift and keeps the other as a
  shift. Candidates on a moved terminus or another N-terminal form than
  the result's cannot be expressed as a mass shift list and are left out.
  PrSM XML: `<alternative_mass_shift_lists>` with one `<mass_shift_list>`
  per alternative under the proteoform (read back by `parseXml`; older
  files have none). PrSM TSV: the column "Alternative PTM explanation"
  after MIScore (`getAlternativePtmStr`: each alternative's variable PTMs
  and unexpected shifts in position order, `Name:[a-b]` / `+mass:[a-b]`
  joined by `;`, the alternatives joined by ` | `). The SQLite PrSM table
  does not have it.
- **Research aids.** `TOPPIC_LOCAL_ALIGN_PARA` ("name=value,..."): `iso`
  (0; 1 = the isotope match of `stat/local`), `max_ptm` (2; 1 = single
  PTMs only), `split` (1), `term` (1), `term_range` (5),
  `term_max_matched` (3), `tie` (1; 0 = no tie-break rules).
  `TOPPIC_LOCAL_ALIGN_LOG=<file>`: one tab-separated line per shift
  handled (spectrum id, PrSM id, protein, shift index and mass, piece with
  the break site and the matched sites of the range, the explanations,
  the decision small / none / ambiguous / localized / no_sites / overlap,
  the sites with the confidence and the score).
- **Results** (toppic_test experiment 140, report 0140; 480_F3_02 in
  setting d plus `-B` with Oxidation M, Acetyl K, Methyl KR, Dimethyl KR,
  Trimethyl K, Phospho STY; the same 3,354 / 1,997 / 438 identifications
  with equal E-values in every run, the stage 2-3 s): the branch replaces
  400 shifts of 311 PrSMs (main 278 of 265): 222 one-shift PrSMs (116
  Phospho, 42 Phospho+Phospho, 36 Methyl ...) and 109 pieces of the 163
  two-shift PrSMs (76 Oxidation; 20 PrSMs with both shifts explained, 69
  with one); 44 one-shift PrSMs are ambiguous (Dimethyl / Methyl+Methyl
  19, Acetyl / Trimethyl 16), 42 localized shifts match one isotope away,
  29 PTMs are on the fallback range. 148 PrSMs are localized by both
  methods, 145 with the same PTM(s) and 141 with the same best site(s).
  Main alone localizes 117 (76 by truncated candidate forms, 12 by
  dropping the N-terminal acetylation form, 29 by choosing among
  ambiguous explanations); the branch alone 163 (63 one-shift PrSMs by
  pairs, 24 by single PTMs that main's filters rejected, 86 pieces of
  two-shift PrSMs).

## toplib_build (`src/toplib`, `console/toplib_build`)

toplib_build is a C++ port of TopLib's `src/ms_library_building.py`
(`~/claude/toplib`, commit 7577b09): it clusters identified MS/MS spectra
and writes target and decoy representatives into the SQLite table
`target_decoy_spectra_representatives`. With the hidden
`--python-compatible` the port reproduces the Python script's library value
for value, and the tests that show it (toppic_test experiments 81 and 82,
reports 0081 and 0082: data/480_lib and 480_F3_02's TopPIC 1.9 and 2.0
output) compare every cell. By default it departs from the script in two
ways, and `--disable-different-proteoform-merge` adds a third (experiment
82):

- **A proteoform is the pair (proteoform id, protein accession)** in the
  step that keeps one cluster's identification per proteoform and
  precursor charge: TopPIC 1.9 gave PrSMs of different proteins sharing a
  precursor feature one proteoform id (43 of 1,818 ids on 480_F3_02; 2.0's
  cluster protein correction removed that), and the script, which groups
  by the id alone, dropped one protein's identification when two proteins
  shared an id.
- **The PrSM table's optional columns are also found under TopPIC 1.9's and
  2.0's names** (`idColumnMap`: TopLib's name, TopPIC 1.7's, first; then
  e.g. `Unexpected modifications`, `# Masses`, `# Matched masses` /
  `# Matched exp fragments`). The script looks for 1.7's names only, so
  from 1.9 and 2.0 tables it copied none of the eight renamed columns.
- **`--disable-different-proteoform-merge`** (visible, argument key
  `disableDifferentProteoformMerge`, default false; `build`'s
  `disable_different_proteoform_merge`): spectra of different proteoforms
  (the same pair) are never clustered together. In the spectral
  clustering a pair of identified spectra of different proteoforms gets
  the distance n * n (n the precursor cluster's size) instead of the
  cosine distance (`clusterGroup`'s `cannot_link`): every average linkage
  merge containing such a pair then averages at least 4, above the 0.7
  cutoff, so no flat cluster holds it, while unidentified spectra can still
  join any cluster. This was the default up to main 47b9fc57; the user made
  the script's behaviour the default (2026-10-04). By default every pair
  gets its cosine distance, so the clusters are the script's: they join
  different proteoforms (46 of 47 multi-member clusters on data/480_lib,
  34 of them different proteins), and as each representative keeps one
  identification, 50 of 4,630 proteoforms get none (8 of 1,908 on
  480_F3_02's TopPIC 2.0 output, 21 of 1,868 on its 1.9 output). With the
  option the library and cluster files equal main 47b9fc57's default byte
  for byte on experiment 82's inputs; without it the cluster files equal
  `--python-compatible`'s. The parameter block prints
  `Different proteoform merge: enabled|disabled`. Do not replace the
  option by a stricter cosine cutoff: that changes the clustering of every
  spectrum, not only of the mixed ones, and still does not guarantee the
  separation.

`--python-compatible` turns the first two off (it is what the comparisons
with the Python script use); it is independent of
`--disable-different-proteoform-merge`, so without that option it
reproduces the script. Facts to keep straight:

- **The Python choices are kept on purpose**, including odd ones: the
  top-50 peaks, precursor clusters by charge then complete linkage within
  2.2 Da, the 10 ppm cosine distance with an isotope shift (TopPIC's
  `mass_constant::getIsotopeMass()`, 1.00228 Da; the script's 1.00235 Da
  is not kept) and its
  early return when the sum passes 1, average linkage at 0.7, the
  `weighted_mass_sum_inte_cal` merge (a peak of the first spectrum takes
  every matching peak of the second; a second-spectrum peak may be merged
  into several), the precursor representative from the truncated integer
  mass differences, the lowest-E-value identification and the removal of
  all but the best cluster per (precursor charge, proteoform id), the decoy
  loops (redraw while negative, then while above the precursor mass, the
  second loop not rechecking the first), and pandas' column types. Each
  function's comment names the Python function it ports; keep that.
- **Exact arithmetic of the Python libraries** (`lib_numeric`,
  `lib_cluster`): SciPy's `nn_chain` linkage (with its tie rule, the stable
  sort by distance and `label`), `fcluster(criterion='distance')`'s
  traversal order (the cluster ids), numpy's pairwise sum, CPython's
  `sum()` (Neumaier, Python >= 3.12) for the PPIR, CPython's Mersenne
  Twister with `random.seed(int)` and `random.sample` (the decoys), pandas'
  `precise_xstrtod` for the PrSM table's numbers (not always correctly
  rounded: 395 of the data set's 5,360 E-values differ from `strtod` by an
  ulp), and `np.linalg.norm` as OpenBLAS's AVX-512 (SkylakeX) `ddot`
  computes it, with `std::fma`. The `src/toplib` sources are compiled with
  `-ffp-contract=off` (CMakeLists.txt) so that nothing else is fused.
- **Not reproduced:** numpy's `log2` (Intel SVML on AVX-512) differs from
  glibc's by an ulp for about 0.006% of values, so a few normalized
  intensities differ by an ulp; numpy's default sort is not stable and the
  port uses stable sorts, so equal values (two peaks of one mass with
  different charges, equal E-values) may come out in another order. The
  norm matches numpy only where numpy's OpenBLAS picks the SkylakeX kernel.
  Everything else was bit-identical in experiment 81.
- **Inputs:** a TopFD msalign file plus a PrSM table (the script's input;
  the header is the first line with a `Spectrum ID` column, so TopPIC's
  parameter block is skipped, which TopLib does with a separate script), or
  a library msalign file (v1.2.1 keys, identifications in the headers;
  `lib_msalign_reader::readLibraryMsalign` maps them onto the table's
  columns). The optional table columns are read under TopLib's names or
  TopPIC 1.9's and 2.0's (see above).
- The hidden option `--cluster-file <file>` writes every spectrum's cluster
  (for comparing the clustering with the Python pipeline's, which it does
  not write).

## toplib_search (`src/toplib/lib_search*`, `console/toplib_search`)

toplib_search is a C++ port of TopLib's `src/ms_library_query.py` (the
same commit): it searches the spectra of a TopFD msalign file against the
library table and writes `toplib_output/query_res.tsv`. It reproduces the
script's table byte for byte (toppic_test experiment 83, report 0083:
480_F3_02 against five libraries of data/480_lib and 480_F3_02, the
parameter variants and a stress query); there is no compatibility switch
because nothing in the search was changed. Facts to keep straight:

- **The script's choices are kept on purpose**: the library rows with a
  proteoform, targets and decoys alike; only the query's first precursor;
  the ppm test in ppm of the *library* mass with the isotope shift
  (`mass_constant::getIsotopeMass()`, 1.00228 Da, not the script's
  1.00235; Da mode has no shift); the cosine loop of
  `peak_50_cosine_mass_ppm` with the library spectrum first and its early
  return; the *first* candidate of the smallest distance (`np.argmin`,
  a NaN distance first) kept when Python's `min()` is at most 1; the
  cutoffs as `np.arange(0.1, 1.0, 0.001)` computes them (start + i *
  (the second value - start), not i * 0.001); the cutoff whose FDR is
  *closest* to 1% (it may be above 1%); the similarity filter of 0.3
  afterwards. The decoys share their target's precursor mass and rarely
  win (40 of 3,036 best matches on 480_F3_02), so the chosen distance
  cutoff is about 0.98 and the 0.3 similarity filter, not the FDR,
  decides which matches are written (none of the 40 decoys reaches it).
- **pandas' output**: the library columns keep the dtype `read_sql_query`
  gives them over the whole table (integers only -> int64; any REAL value,
  or an integer column with a NULL -> float64, written `924.0`; text or
  all NULL -> object), floats are written as Python's repr
  (`lib_numeric::pythonRepr`, `std::to_chars` shortest digits laid out as
  CPython's `format_float_short`), NULL as an empty cell, text quoted as
  `csv.QUOTE_MINIMAL` quotes.
- **Deliberate differences**, none of which changes a result on a single
  TopFD file: the library path is used as given (the script prepends
  `./TopLib/`); `-T` and `-c` are case-insensitive and checked (the
  script's `-T da` takes the 2.2 default but tests in ppm); tolerances
  must be nonnegative numbers; spectra of a repeated `SPECTRUM_ID` (a
  combined msalign file) are searched and reported one by one, where the
  script merges their masses and then fails ("All arrays must be of the
  same length").
- **Not reproduced**, as in toplib_build: numpy's `log2` (SVML) differs
  from glibc's by an ulp for a few values, which moves a distance by an
  ulp (1 query of 4,257 on 480_F3_02); pandas' descending `sort_values`
  and `np.argsort` are not stable, so a tie of intensities across the
  50-peak cut or two equal masses may be ordered differently (stable sorts
  here; none in 480_F3_02). The norm is `lib_numeric::numpyNorm` (OpenBLAS
  SkylakeX `ddot`), and `-ffp-contract=off` keeps the loops unfused as
  numba's are.
- The library search needs no index: candidates come from the library
  sorted by precursor mass with a window that covers the exact test, then
  in library order; a search of 4,257 spectra takes about 0.5 s.
- The hidden option `--match-file <file>` writes every query's match
  before the FDR step (`spectrum_id`, `file_name`, `lib_index`,
  `representative_id`, `flag`, `calculated_cosine_dist`), the table the
  script computes in `ms_query` and does not store.

## Isotope distributions (`common/emass`)

`common/emass` is a C++ port of toppic_test's `emass.py` (experiments 107,
108 and 111): emass's convolutions of isotope patterns (abundance and
abundance-weighted mass offset per number of extra neutrons, trailing peaks
<= 1e-100 pruned, binary powers for whole atoms) and J.C.P. Miller's
recurrence for fractional atom numbers (`BelowOne::kMixture` takes fewer
than one atom as a mixture with and without it). `isotope_data` holds the
isotopes of C, H, N, O and S (AME2020 masses, IUPAC representative
abundances) and the averagine compositions: `kHumanAveragine` (the UniProt
human reference proteome, experiment 108; the default),
`kHumanNoSulphurAveragine` (the same proteome without cysteine and
methionine, so C, H, N and O only; experiment 112) and `kSenkoAveragine`. `averagine_table::writeAveragineTable` writes TopFD's
`res/base_data/theo_patt.txt` and reproduces the file of experiment 111
byte for byte (toppic_test experiment 112, in 2.5 s against `emass.py`'s
14 s). Facts to keep straight:

- **Bit-identical to `emass.py`**, which experiment 112 checked for every
  rounding and fractional atoms on 14 masses from 100 Da to 100 kDa: the
  monoisotopic masses and the averagine unit mass are summed with
  Neumaier's compensation as Python's `sum()` does (Python >= 3.12), the
  convolutions sum each output in the order of the longer series as
  `numpy.convolve` does, and `averagineFormula` rounds halves to even as
  Python's `round()`. The elements are always convolved in the order C, H,
  N, O, S; `emass.py`'s H-filled formulas list H last, which moves those
  distributions by about 1e-16.
- **It does not replace `EnvBase`.** TopFD still reads the table from
  `theo_patt.txt` (`ms/env/env_base`); `common` cannot depend on `ms`, so
  the table's 11,000 entries and 10 Da interval are repeated in
  `averagine_table.cpp`. Keep the two in step.
- **`theo_patt_gen`** (`src/console/theo_patt_gen.cpp`) is the command-line
  front end: `theo_patt_gen [-o file] [-c human|human_no_sulphur|senko]
  [-s k]` (`--output`, `--composition`, `--sulphur-num`; defaults
  `theo_patt.txt`, `human`, 0). It is a
  developer tool, declared
  `EXCLUDE_FROM_ALL`: `make` does not build it (`make theo_patt_gen`
  does), it is deliberately absent from `TOPPIC_EXECUTABLES` (CMake
  documents installing an `EXCLUDE_FROM_ALL` target as undefined behaviour)
  and it has no manual. Regenerate the resource with
  `bin/theo_patt_gen -o res/base_data/theo_patt.txt`. `--sulphur-num k`
  (0-5; k > 0 only with `human_no_sulphur`) gives every formula exactly k
  S (`averagineFormulaWithSulphur`): the composition is scaled to the entry
  mass minus k 32S masses (at least 0, so the lightest entries are S_k
  alone) and k S are added.

## Index files version 2, the combined filter pass and the E-value stage (toppic_test experiment 150)

Branch `speedup_four_ptm_s2` (toppic_test experiment 150, report 0150;
from main 915bb387, merged into main on 2026-10-09). `toppic -u 16 -f C57
-d -t FDR -T FDR -y FDR -s 2 -b topmg_var_mods.txt` on 480_F3_02 takes
543 s instead of 816 s (-33%) with every output file byte-identical to
main 915bb387's (the `-k` stage files, the record PrSM file, the PrSM
XML). Facts to keep straight:

- **Index files version 2** (`MassMatch::serializeMassMatchV2`): the
  header `TPMMIDX2` (scale, proteoform / column / row numbers), then the
  seven arrays as zstd -1 frames, the column bounds stored once
  (`col_index_bounds_`, col_num + 1 ints) instead of the begins and ends.
  `deserializeMassMatch` probes the magic and reads a version 1 Boost
  archive or a version 2 file (one `read` into a thread-local buffer,
  `posix_fadvise DONTNEED` after it). `topindex --index-format <1|2>`
  (hidden, default 2). The human reviewed database with decoys takes 93 GB
  instead of 172 GB. Version 1 index directories stay readable; the file
  names (`filter/mng/index_file_name`) are unchanged, so the format is
  known only by reading a file. `MassMatch::compScores` prefetches the
  column bounds and row ids a few masses ahead (no measurable change).
- **One pass for the zero-shift, variable PTM and one-shift filters**
  (`filter/varptm/zero_var_ptm_filter_processor`): with a variable PTM
  file and no approximate spectra, `toppic_process` runs one stage, "Zero
  unexpected shift, variable PTM and one unexpected shift filtering" ("Zero
  unexpected shift and variable PTM filtering" with `-s 0`), that loads
  each block's four zero-shift index files once and scores them with the
  `ZeroPtmFilter`, a `VarPtmFilter` built on the loaded indexes (a new
  constructor) and, with `-s` of at least 1, the `OnePtmFilter`. The
  one-shift filter's index files are the zero-shift ones cut at its column
  limit (`truncateColumns`), so it scores the loaded indexes with its
  windows cut there (`MassMatch::compScores(..., col_limit)`), and where a
  spectrum's mass lists equal the variable PTM filter's (both
  `getIntMassErrorList` at precursor offset 0) and every window lies below
  the limit, it finds its candidates on the variable PTM filter's score
  arrays in place (`VarPtmFilter::computeBestMatch`'s `base_scores_hook`,
  called before the precursor bonus is added). The three stages' output
  files are unchanged (`toppic_zero_filter_*`, `toppic_var_filter_*`,
  `toppic_one_filter_*`), so the searches and the tag filter's hooks read
  them as before; the one-shift index files are not read in this mode
  (topindex still writes them for the separate stages and topmg). Hidden
  `--separate-var-ptm-filter` runs the separate stages; without a variable
  PTM file the zero-shift stage runs alone as before. The zero-shift
  filter computes its N-terminal index scores once per spectrum and copies
  them for the three precursor offsets (the N-terminal mass lists do not
  depend on the precursor); the stage's three spectrum readings run in
  parallel threads. Copying the score arrays between the filters was
  tried and cost more than it saved (cache eviction); read them in place.
- **The multiple-shift search's diagonal finder**
  (`CompShiftLowMem::findBestShift`, `findBestShiftChunked`) counts the
  mass differences in chunks of 65,536 bins from the top of the range
  down, keeps the count events of at least a threshold (4, then 3, then 2)
  in the original loop's order per spectrum mass (a window across a
  chunk's lower edge waits in a pending buffer until the lower chunk has
  emitted its part) and replays them through `checkShift`; when the
  replay does not fill the 20-shift list the next threshold is tried, and
  last the original loop runs. The lists are the original's (an event
  below the final list's smallest count never changes it; checked on a
  dump of 5,046 calls, `TOPPIC_SHIFT_DUMP=<file>`); the search went from
  157 to 35 s, because the original loop touched a random cache line per
  mass pair in a buffer of up to 770 MB.
- **The E-value stages** (`stat/tdgf`): a PrSM with at most
  `comp_evalue_min_match_theo_frag_num_` (4) matched fragments gets the
  maximum E-value without a table and is left out of its spectrum's
  groups; a spectrum whose PrSMs are all such skips the generating
  function. `CountTestNum` is built in a background thread started before
  the filters (`EValueProcessor::prepareCountTestNum`, joined when the
  stage starts). The PrSMs are written in spectrum order through the
  `OrderedTaskRunner` (16 pending tasks per thread; fewer let one slow
  spectrum idle the others) with `PrsmXmlWriter::writeTrimmedLines`, the
  bytes the former per-thread files and merge produced.
- **The generating function in score slices** (`CompProbValue`): the
  table (height = the group's best score + 1, at most 128, x the mass in
  units of 1/27.4 Da) is swept in slices of 8 rows (`kChunkRows`), each a
  sweep over all columns in blocks of 1,350 columns through a circular
  page of 6 blocks, which needs the slice below only at the columns with
  a position score (`saved_rows_`); the look-back window of the residues
  stays in L2, where the whole-height sweep (a page of 5,000 Da x 128
  rows, 138 MB per thread) streamed the table through memory and the stage
  was bound by memory bandwidth. The per-cell operations and their order
  are the old sweep's, but every slice of every layer starts from a
  cleared page, so the columns before column 0 hold zeros (branch
  `evalue_clean_sweep`, toppic_test experiment 152). **The old sweep
  carried a layer's last pass into the next layer**: it cleared the tail
  of its page once per computation, so for precursors whose table wraps
  the page (above about 4.6 kDa) a shift layer's first blocks read,
  through the wrap, the last 5,105 columns (the largest residue mass) the
  previous layer left there: random sequences with layer 0's scores at
  some mass near the end of the table, placed at mass 0 of layer 1, which
  then collect the spectrum's matches a second time. Only the one-shift
  and two-shift probabilities of such precursors were raised (the
  zero-shift layer and the variable PTM PrSMs are layer 0): by a median
  of 10^0.00002 at 5-9 matched fragments and 10^0.0001 at 10-19, but by
  up to 10^8.7 for PrSMs with 20 or more (one-shift targets of the first
  pass above 4.6 kDa); decoys by at most 10^0.02. Main 3b8041dd
  reproduced it (`tail_capture_`) to keep main 915bb387's output; on
  480_F3_02 (setting d with the four variable PTMs) the clean sweep gives
  the same 3,412 PrSMs / 2,016 proteoforms / 447 proteins and cutoffs,
  and 5 spectra report another PrSM of the same protein (a one-shift
  explanation that now beats a zero-shift PrSM at a +-1 isotope offset
  or an isobaric variable PTM, in 4 of them with more matched fragments).
- **Research aids**: `TOPPIC_PROFILE=<file>` (`common/util/
  stage_profile`: thread-seconds per named part and stage; the output
  files do not change) and `TOPPIC_SHIFT_DUMP=<file>` (the diagonal
  finder's inputs).
- **Results** (480_F3_02, toppic_test report 0150, one run at a time):
  main 820 / 812 s -> 542 / 543 s, user time 10,880 -> 7,510 s, system
  time 529 -> 128 s, peak memory 18.9 -> 18.4 GB; 3,412 PrSMs, 2,016
  proteoforms, 447 proteins at 1% FDR in every run. Not kept:
  binary-search shift predecessors in the alignments (the quadratic loop
  was 1% of the search), per-class table heights in the E-value.

## Rescoring of the PrSMs (`stat/rescore`, default with `-d`)

Branch `prsm_rescore` (toppic_test experiment 151, report 0151 and its
amendment 1; from main 3b8041dd, merged into main on 2026-10-09 as the
hidden `--rescore`, off by default). **On by default in a target-decoy
search (`-d`) since 2026-10-09**: argument key `rescore`, "true"; the hidden
`--disable-rescore` sets "false" (`--rescore` is gone), and
`validateArguments` sets it to "false" without `-d` and with
`--fdr-groups` (both were errors with `--rescore`), so a run without `-d`
is unchanged; the parameter block prints "FDR computation: rescored PrSMs
(linear SVM, 3-fold cross-validation)" only when it runs (this line and
`--fdr-groups`' "FDR computation: separate PrSM groups" are printed with or
without a variable PTM file; before, only with `-b`). It gives the best
PrSM of every spectrum a score of a semi-supervised linear model, and the
FDRs rank by it instead of the E-value. A `-d` run's FDRs, cutoffs and
tables therefore differ from those of main before the change;
`--disable-rescore` reproduces them. Facts to keep straight:

- **Where it runs.** `TopPIC_postprocess` keeps its order (clustering,
  cluster protein correction, recount, post mass matching, PTM
  characterization) and runs `prsm_rescore::process` on the result
  (`toppic_rescore`) just before the FDR; `TopPIC_multi_postprocess` (`-c`)
  does the same on the merged PrSMs and spectra (not tested). The E-values
  are not changed; `prsm_fdr::process(..., use_rescore)` sorts by
  `PrsmStr::cmpRescoreDecProtInc` and its three FDR functions compare the
  key `rescoreKey` (the score negated; `eValueKey` by default), so the
  cutoffs (`-t FDR`, `-T FDR`, `-y FDR`) apply to the rescored FDRs.
  `-t EVALUE` cutoffs still compare E-values.
- **The score** (`stat/rescore/rescore_model`): Percolator's algorithm
  (Kall et al., Nat. Methods 2007). The spectra are split into 3 folds at
  random (`std::mt19937` with seed 1 on the spectrum ids in increasing
  order, `rng() % 3`); for each fold a linear SVM
  (`stat/rescore/linear_svm`: L2 penalty, squared hinge loss, C = 0.1,
  "balanced" class weights, dual coordinate descent, the bias as a
  penalized feature; scikit-learn's `LinearSVC` defaults) is trained on the
  other two folds 5 times, each time with the targets at 1% FDR of the
  current score (TopPIC's ratio, `passFdr`) as the positives and every
  decoy as a negative, the first iteration ranking by -log10 E-value; the
  features are standardized on the training folds. The held-out fold's
  scores are put on one scale (the training folds' 1% threshold -> 0, their
  median decoy -> -1). A fold with fewer than 10 positives gives every PrSM
  of the search its -log10 E-value (the FDR is then the E-value's). The
  searches (`<tool>` toppic, topmg) are rescored apart, as the FDR counts
  them apart. The `<rescore>` element of the PrSM XML (after
  `<protein_fdr>`, written and read only when set: `Prsm::getRescore`,
  `PrsmStr::getRescore`, NaN without) carries the score to the FDR and
  the outputs; the TSV tables do not show it.
- **The features** (`prsm_rescore::compFeatures`, names in
  `getFeatureNames`): -log10 E-value (clipped to [-5, 100]), log10 test
  number, -log10 one-protein probability, matched fragments and masses (the
  record's, after the recount and the post mass matching), matched masses
  per spectrum mass and fragments per residue, log10 masses of the
  spectrum, 0 or 2 unexpected shifts, variable PTMs, log10(1 + total
  |shift|), a shift within 0.02 Da of a common modification mass, NME,
  acetylated N-terminal form, start <= residue 1, end at the protein's
  last residue, log10 length and mass, precursor charge, log10 intensity,
  intensity / the largest precursor's, precursor number, log10 feature
  intensity, ECScore; and from the b / y (c / z) pairs of
  `peak_ion_pair_util::genePeakIonPairs` on the spectra read with the
  EnvCNN cutoff of `-F` (the masses the search used): matched intensity
  fraction, median and median absolute deviation of the ppm errors (pairs
  at an isotope offset left out), |median|, isotope-offset share,
  cleavage-site coverage, the longest run of covered sites, complementary
  sites. One more: |RT - predicted RT|, the prediction a ridge regression
  (alpha 1) on the amino acid composition, log10 length and log10 mass,
  fitted in each fold on the training folds' targets at 1% FDR by E-value.
  The console prints, per search, the PrSMs, decoys and the targets at 1%
  FDR by E-value and rescored, and the 8 largest standardized weights.
- **The FDR caveat** (toppic_test report 0151's follow-up): with this
  order the cluster protein correction has already given some spectra the
  proteoform of their cluster's best PrSM (on 480_F3_02 331 records, 140 of
  them decoys before the correction, 192 at E = 1e300), so the decoy count
  does not check those records, and a rescored PrSM among them is accepted
  on the cluster's identification (feature-level propagation) rather than
  on its own spectrum. Rescoring the PrSMs before the correction keeps the
  decoys valid (+2.4% PrSMs instead of +6.3% in the Python prototype).
- **Two rounds of FDRs** (branch `rescore_fdr_fallback`, toppic_test
  experiment 153, report 0153; from main 2e97163f, not merged).
  `prsm_fdr::process` with `use_rescore` computes the three FDRs by the
  E-value (on copies of the PrSMs, sorted from the input order as the
  E-value computation sorts them, so the E-value round equals a
  `--disable-rescore` computation) and by the rescore (on the PrSMs, as
  before), counts the targets within each level's cutoff (`FdrCutoffs`:
  the `-t` / `-T` / `-y` value when the type is FDR, else 1%;
  `getFdrCutoffs` in `toppic_process`) and keeps a level's rescored FDRs
  unless the E-value's give more targets (a tie keeps the rescore; each
  search chooses on its own; the console prints "FDR computation -
  <tool>: targets within the cutoffs by E-value / rescore: ..."). Only the
  PrSMs whose E-value FDR differs are rewritten, so where the rescore wins
  at every level the files are those of a single rescored round. The
  hidden `--disable-rescore-fallback` (argument `rescoreFallback`) keeps
  the rescored FDRs at every level; the parameter block's FDR line ends
  with ", the E-value at a level where it gives more" when the rule is on.
  It is used by `TopPIC_postprocess`, `TopPIC_multi_postprocess` and the
  second precursor merge. Results: 480_F3_02 setting d keeps the rescore at
  every level (3,573 / 2,041 / 447 against 3,374 / 2,004 / 440 by E-value)
  and writes main's record, XML and table files but the parameter line;
  the Bruker run of experiment 152 (b_p1) 1,276 / 82 / 71 (PrSMs by the
  rescore, proteoforms and proteins by the E-value; main 1,276 / 64 / 59,
  `--disable-rescore` 806 / 83 / 72); the 500-spectrum subset of 480 with
  `--mass-calib` (6 decoys) 474 / 345 / 152 by the E-value at every level
  (main 26 / 22 / 12). With `--disable-rescore` the files are main's.
- **The rescore's PrSM-level FDR is not valid on the Bruker run**
  (toppic_test experiment 153, entrapment search: human SwissProt plus 0.5
  x shuffled `ENTRAP_` proteins): at 1% FDR the rescore accepts 1,357
  PrSMs with 59 entrapment PrSMs (estimated FDP 13.0%), the E-value 777
  with 3 (1.2%); 47 of the 59 have E = 1e300 (records the cluster protein
  correction gave a cluster's proteoform and that then failed the E-value),
  and so do 465 of the 1,357. The model clips -log10 E at -5, so such a
  record is ranked by its other features (feature intensity and the like,
  shared with its cluster); a decoy of this kind (rescore 1.06) is also
  what stops the rescored proteoform list after 64 targets. The two rounds
  fix the proteoform and protein levels but cannot see this at the PrSM
  level, where the rescore gives more targets. On 480_F3_02 the rescored
  FDP was 0.84% (report 0151).
- **The lowest rescore for the PrSMs with the maximum E-value** (branch
  `rescore_evalue_floor`, toppic_test experiment 155, report 0155; from
  `rescore_fdr_fallback` e10d3fd0, not merged). `prsm_rescore::process`
  takes an `EValueFloor` (`getRescoreEValueFloor` in `toppic_process`, for
  the three callers): the PrSMs with E >= `ExpectedValue::getMaxDouble()`
  (1e300) are left out of the training and get the lowest double,
  `std::numeric_limits<double>::lowest()` (-DBL_MAX; the user's choice,
  amendment 1, fc2e8b46; f3d913c6 gave 1 below the lowest rescore of the
  search's other PrSMs, the same order), so the FDRs rank them last, as the
  E-value does; a search without a model keeps the -log10 E-value. Every
  use of a rescore is a comparison or a negation (`rescoreKey`), so the
  value is safe; `Prsm::appendXml` writes a rescore of magnitude 1e100 or
  more with `max_digits10` digits (`-1.7976931348623157e+308`, read back
  exactly by `std::stod`), not with `str_util::toString`, whose fixed
  notation writes 309 digits. The console prints "<n> PrSMs with the
  maximum E-value (<d> decoys) get the rescore -DBL_MAX (the lowest double)
  and are left out of the training", and the parameter
  block's FDR line "(linear SVM, 3-fold cross-validation; the maximum
  E-value last)". Hidden `--disable-rescore-evalue-floor` (argument
  `rescoreEvalueFloor`) turns it off (the files of `rescore_fdr_fallback`,
  checked on the Bruker run and the 500-spectrum subset); research
  switches in `TOPPIC_RESCORE_PARA`: `floor_train=1` (the model trained on
  every PrSM as before, only their scores replaced) and `floor_value=<v>`.
  With -DBL_MAX every output equals f3d913c6's but the `<rescore>` values of
  those PrSMs (toppic_test experiment 155, amendment 1).
  Not 0: the rescore's scale puts the training folds' 1% threshold at 0 and
  the median decoy at -1, so a 0 sits at the acceptance threshold and,
  with the tie rule of the FDR (tied decoys first), stops every list there
  (the Bruker run's 2% and 5% lists equal its 1% list). Results at 1% FDR,
  entrapment FDP (human + 0.5 x shuffled `ENTRAP_`) in brackets: 480_F3_02
  setting d 3,481 / 2,046 / 443 (main's rescore 3,573 / 2,041 / 447,
  E-value 3,374 / 2,004 / 440); on the entrapment database 3,475 PrSMs
  (0.9%), 2,032 proteoforms (0.9%), 442 proteins (0), E-value 3,353 (1.0%) /
  1,981 (1.1%) / 439 (2.7%), main's rescore 3,556 (0.8%) / 2,020 (0.4%) / 442
  (0): valid, gain +3.6% PrSMs, all gained PrSMs with E > 1 (the E = 1e300
  ones gone). The Bruker run: 905 PrSMs (rescore) / 82 / 71 (E-value); on
  the entrapment database the floored rescore wins at every level, 918
  PrSMs (3.9%, against 13.0% without the floor and the E-value's 1.2%), 92
  proteoforms with 4 entrapment ones, 80 proteins with 4: still not valid.
  Its accepted entrapment PrSMs have E 0.4-346 and sit on intense features
  (the model's largest weight is the feature intensity once the E = 1e300
  records are out of the training): among the records with E < 1e300, junk
  on features above 1e7 is entrapment and decoy 10 : 10, elsewhere about
  1 : 4 (the expected 1 : 3). Trained on every PrSM (`floor_train=1`) 903
  PrSMs (5.3%), 87 proteoforms with 1; with 0 (`floor_value=0`) 895 (4.4%).
  Few decoys: st_1 (63 decoys) main's rescore keeps 16 PrSMs against the
  E-value's 177, the floored rescore 64, the two rounds the E-value's 177 /
  56 / 35; the 500-spectrum subset of 480 the E-value at every level
  (with `--mass-calib` 481 / 343 / 151, the PrSMs by the rescore). No time
  cost.
- **Results on 480_F3_02** (setting d of toppic_test report 0062, TopFD
  main 915bb387 output of experiment 151; toppic_test experiment 151,
  amendment 1): 3,573 PrSMs / 2,041 proteoforms / 447 proteins at 1% FDR
  against main 3b8041dd's 3,374 / 2,004 / 440 (the toppic_test Python
  prototype, `rescore.py` with the same features: 3,588 +- 4 / 2,045 / 443
  over 5 seeds); 242 spectra gained and 43 lost, 98 of the accepted PrSMs
  with E = 1e300 (cluster-corrected records, the caveat above). The step
  takes a few seconds (12:26 against 12:53 for the run, the same memory).
  Without the option the branch writes main's 62 files (56 equal, 6 but for
  the time lines). On the human + 0.5x shuffled entrapment
  database of experiment 151 (`ENTRAP_` proteins, indexed by this
  topindex in 3:19, 140 GB) the estimated false discovery proportion
  (entrapment PrSMs x (1 + 1 / 0.5) / accepted PrSMs) is 0.84%, 1.66% and
  2.86% at 1%, 2% and 5% FDR (3,557 PrSMs at 1%), against the E-value's
  0.98%, 1.41% and 3.56% (3,353); with 8-20 entrapment PrSMs the estimate
  has a standard error of about 0.3 percentage points at 1%.

## The second precursor search (`prsm/prsm_second_prec`, `--second-precursor`)

Branches `second_precursor` (toppic_test experiment 151, report 0151 and its
amendment 2; from main d31939d9) and `second_precursor_id` (amendment 3: the
units keyed by spectrum id and precursor id, the `Precursor ID` column; from
main 2863fa08). Off by default; toppic's hidden
`--second-precursor` (argument key `secondPrecursor`, "false"; the parameter
block prints "Second precursor search: yes" only then) searches the spectra
with two or more precursors a second time with the second precursor, and
reports up to one PrSM per precursor of a spectrum. Refused with `-c` and
with `-r` other than 1. Facts to keep straight:

- **Why.** TopFD writes every MS1 feature it assigns to an MS/MS spectrum's
  isolation window as a precursor (`PRECURSOR_MZ=a:b:c`, by intensity;
  `deconv_ms2_process.cpp` step 3; 2,440 of the 4,323 spectra of
  480_F3_02 have two or more), and the search uses the first
  (`MsHeader::getFirstPrec*`, `spectrum_set_factory`). A co-isolated
  proteoform's fragments stay unexplained. TopMPI (Wang, Xiong and Liu,
  Proteomics 2025) searches the second precursor after removing the first
  proteoform's matched masses; this is that search on TopFD's precursor
  list.
- **The passes** (`TopPIC_second_precursor` in `console/toppic_process`,
  after `TopPIC_postprocess` and before `TopPIC_report`; `TopPICProgress`'s
  identification stages are `TopPIC_identifyAll`, which both passes call).
  The first pass is a normal run up to its record file. Then
  `prsm_second_prec::writeResidualSpectra` writes `<base>_resid_ms2.msalign`
  from the spectrum file the first pass searched: the spectra with two or
  more precursors, their precursor list with the second precursor first
  (the ids of the scan's list and the feature ids kept: `MsAlignWriter`
  adds a `PRECURSOR_ID=1:0:...` line to a spectrum whose precursors are not
  at the positions of their ids, `MsAlignReader` takes the ids from it, and
  without the line numbers the precursors by position, so TopFD's files
  and every other msalign file are unchanged; the second pass's PrSMs thus
  carry `precursor_id` 1), without the masses that the
  spectrum's accepted first-pass PrSM (a target under the `-t` cutoff)
  explains: within the fragment tolerance of one of the proteoform's
  fragment masses (`theo_peak_util::geneProteoformTheoPeak`, so the shifts'
  ranges give no ions and the N-terminal label counts), as it is or one
  isotope (`mass_constant::getIsotopeMass`) lower or higher. The feature
  file is copied to `<base>_resid_ms2.feature`. The second pass runs the
  same identification stages and post-processing on that file with its own
  copy of the arguments (`rescore` "false": the merged set is rescored
  instead), so its files are `<base>_resid_ms2.toppic_*` and, post-matched,
  `<base>_resid_post_ms2.*`: `prsm_post_mass_match::sqlFileName` strips
  `_resid` (and then `_calib` and `_mono`) to find `<base>.sqlite`, and
  `postBaseName`
  keeps `_resid` in the post-matched name.
- **The units.** A PrSM is the unit (spectrum id, precursor id), the
  precursor id the position of its precursor in the scan's list (0 the
  first pass's, 1 the second pass's); the spectrum ids are the scan's own
  (until 2026-10-09 the second pass's spectra were appended to the
  spectrum file with their ids raised by `SpPara::getMaxSpecNumPerFile()`,
  10,000,000, and a raised id marked a second unit). The first pass's
  record is kept as `toppic_first_record_prsm`, and `mergePrsms` writes
  the units to `toppic_second_merged`: every first-pass PrSM, and every
  second-pass PrSM whose protein (decoy prefix aside) and precursor feature
  no first-pass PrSM of the spectrum has; one that repeats a first-pass
  PrSM's protein or feature is kept in its place only with a better
  E-value (a truncated or modified form of the first proteoform is not a
  second unit: the residual fragments can still be its); in spectrum id,
  then precursor id order. The units are clustered as `TopPIC_postprocess`
  clusters (`prsm_feature_cluster`, or `prsm_simple_cluster` with `-x`;
  `toppic_second_cluster`, `toppic_second_prot`), without a cluster
  protein correction, recount or post mass matching (each pass did its
  own), rescored (default with `-d`, not with `--disable-rescore`;
  `toppic_second_rescore`) and given one FDR over all units
  (`toppic_second_fdr`; `--fdr-groups` applies), which becomes the file's
  `toppic_record_prsm`. `TopPIC_report` then runs as usual.
- **The spectra of a unit.** `prsm_reader_util::readPrsmsWithSpectra`
  attaches the spectra of a PrSM's spectrum id, every PrSM of the id from
  the one spectrum file; with its `second_sp_file_name` (the residual
  post-matched spectra, `<base>_resid_post_ms2.msalign`, whose first
  precursor is the scan's second) a PrSM gets, from either file, the
  spectrum whose first precursor has its precursor id, as
  `prsm_util::isMatchMs` matches. The rescoring of the merged set reads
  with the residual file, so a second unit's fragment features are
  computed on the masses its search used, as its matched numbers were; the
  report's tables read the one output file, so a second unit is shown on
  the scan: `# Masses` is the scan's, the matched numbers the search's,
  and `Charge` (`MsHeader::getPrecPtrById`, the first precursor when none
  has the id) is the unit's precursor's. The output spectrum file is the
  first pass's post-matched file with the masses the second pass's post
  mass matching added (`mergePostMatchedMasses`: a residual mass the scan
  does not have within the strict fragment tolerance, appended with the
  next peak ids; not with `-D`), one spectrum per scan with the scan's
  precursor list, so its spectrum ids are TopFD's and the SQLite `prsm`
  rows join `ms2_spectrum`.
- **The outputs name the unit.** `Precursor ID` after `Spectrum ID` in the
  PrSM, proteoform and protein tables of `PrsmMatchTableWriter` (toppic
  and topmg; 0 in every run without the option; the rows of a PrSM's other
  protein matches have an empty field), `precursor_id` after `spectrum_id`
  in the SQLite tables `prsm`, `proteoform` and `protein`, and
  `<precursor_id>` in the PrSM XML, which it always had.
- **Cleanup.** `cleanToppicDir(..., copy_outputs = false)` deletes the
  residual pass's intermediates (no XML outputs to copy), and without `-k`
  the residual spectrum, feature and index files too; the merge's files are
  in its lists. The residual pass's record file stays, as record files do.
- **Not supported:** combined runs (`-c`: the merge works per file and the
  fractions' ids would collide), several spectra per set; `--topmg` modes
  are run in both passes but were not tested.
- **The FDR caveat.** The decoys of the merged set are each pass's own, so
  the target-decoy count checks a second unit against the second pass's
  decoys on the same residual spectrum. What it cannot check is a second
  PrSM explained by the first proteoform's fragments; the mass removal is
  what limits that (toppic_test report 0151 section 4.2: 1% of the second
  PrSMs share most of their matched masses with the first against 23%
  without the removal). A spectrum's two units are not independent
  measurements of two proteoforms when the precursors are one species'
  isotopes or adducts; TopFD's features separate those by mass.
- **Results on 480_F3_02** (setting d of toppic_test report 0062, TopFD main
  915bb387 output; toppic_test experiment 151, amendment 2): 2,440 spectra get a second
  precursor; 61,231 masses of 1,934 identified spectra are removed; the
  second pass finds 2,357 PrSMs (targets and decoys), of which 2,103 join
  the first pass's 4,076 units, 40 replace a first-pass PrSM of the same
  protein or feature and 214 are dropped for one. 4,100 PrSMs / 2,382
  proteoforms / 529 proteins at 1% FDR against main d31939d9's 3,374 /
  2,004 / 440 (825 PrSMs of the second pass), the prototype's combination
  of two runs (toppic_test report 0151, `resid_prec.py` + `combine.py`)
  giving 4,107 / 2,394 / 528 with 4,098 accepted units in common; with
  `--rescore` (now the default with `-d`; the first pass rescored too, so the residual spectra follow
  the rescored acceptance) 4,507 / 2,559 / 565 (prototype 4,517 / 2,565 /
  554). 25.5 min against 12.5 (the second pass is a full search of 2,440
  spectra), 7-8.6 GB. Without the option the 62 output files equal main's
  (6 but for the time lines). On the human + 0.5x shuffled entrapment
  database of experiment 151 (one run with both options, 35.7 min, 10.6 GB)
  the estimated false discovery proportion (entrapment PrSMs x (1 + 1 / 0.5)
  / accepted PrSMs) by E-value is 0.96%, 2.01% and 4.92% at 1%, 2% and 5%
  FDR (4,044 PrSMs at 1%, 786 of the second pass, 8 of the 13 entrapment
  PrSMs of the second pass), and rescored 0.40%, 1.26% and 3.72% (4,459),
  the prototype's 0.96% and 0.80%; with 6-28 entrapment PrSMs the estimate
  has a standard error of about 0.3 percentage points at 1%.
- **Results with the units keyed by (spectrum id, precursor id)** (branch
  `second_precursor_id`, toppic_test experiment 151, amendment 3; the same
  TopFD output): the same 4,100 / 2,382 / 529 accepted units with
  `--disable-rescore` (825 of the second pass), rescored 4,519 / 2,563 /
  568 against the raised-id run's 4,507 / 2,559 / 565 (4,493 units in
  common; the folds are drawn over the spectrum ids, which a scan's two
  units now share); 8,274 masses of the second pass's post mass matching
  added to 1,831 scans; 25-26 min, 8.0-8.3 GB. Without the option the run
  writes main 2863fa08's files but the `Precursor ID` column: on the first
  500 spectra 57 of 63 files are identical, the 6 tables equal without the
  column and the time lines, the SQLite tables equal without
  `precursor_id` (0 in every row).

## Two-pass m/z recalibration (`ms/calib`, `topcalib`, `toppic --mass-calib`)

Branch `mass_calib` (toppic_test experiment 152, report 0152 and its
amendments 2 and 3; from main 2863fa08, merged into main on 2026-10-09 as
140d9479: on the 500-spectrum subset of 480_F3_02 in setting d the merged
main writes main 6949cad7's files without the option and the branch's model,
spectra and PrSM files with `--mass-calib`). `topcalib` is the C++
form of the experiment's `calib.py` and `recal_msalign.py`
(`tests/202610/scripts/two_pass_recal/` of toppic_test; the manual is
`docs/manual/topcalib_manual.md`), and toppic's hidden `--mass-calib` runs
the method inside a search. Facts to keep straight:

- **The method.** A zero-shift search (`toppic -s 0`, `-D` is enough) gives
  the calibrants: PrSMs with E-value <= 1e-5, no unexpected shift, no
  variable PTM (`mass_calib_fit::readCalibrants`). Precursor errors
  (precursor mass - proteoform mass, ppm of the proteoform mass) of the
  PrSMs whose difference rounds to 0 isotopes of 1.00235 Da; fragment
  errors of the spectrum's masses within 20 ppm of the proteoform's N- and
  C-terminal ion masses (the activation's ion types, the nearest one). Two
  models fitted apart by robust least squares (`fitModel`, `robustFit`: 10
  rounds of weights sqrt(min(1, 1.5 s / |r|)), s = 1.4826 median |r|,
  weighted normal equations): MS1 e = c0 + c_rt (rt - rt0) / 60 on the
  precursors, MS/MS e = c0 + c_mz (m/z - 800) / 1000 + c_rt (rt - rt0) /
  60 on the fragments; rt in minutes, rt0 the median of all calibrants.
  The two levels are fitted apart because their offsets differ: on the
  Bruker maXis run of experiments 151 and 152 MS1 -3.16 ppm, MS/MS -4.17
  ppm.
- **The shortcut** (`msalign_mass_calib::process`): the correction is
  applied to TopFD's msalign masses, not to the raw peaks: m/z' = m/z (1 -
  e 1e-6), so a neutral mass M of charge z becomes M - e 1e-6 (M + z
  proton). Precursors with the MS1 model at their m/z, fragment masses with
  the MS/MS model at (M + z proton) / z, both at the spectrum's
  `RETENTION_TIME`; six decimals; every other line copied. Experiment 152
  found it equal to recalibrating the raw peaks and rerunning TopFD at
  `-e 5` (shortcut against raw peaks: Bruker 851 / 850 PrSMs, 480 3,308 /
  3,310, the same targets at E <= 1e-5). The SQLite peaks that TopPIC's post mass matching reads are
  not corrected.
- **The model file** (`MassCalibModel::write`/`read`) has calib.py's keys
  (`mz0`, `rt0_min`, `ms1` / `ms2` with `c0`, `c_mz`, `c_rt_per_hour`, the
  calibrant numbers, the robust SDs); it is written by hand with `%.17g`
  (Boost's property tree writes numbers as strings) and read with Boost's
  JSON parser, so either tool reads the other's file. `apply` needs only
  `mz0`, `rt0_min`, `ms1`, `ms2`.
- **Masses.** topcalib takes the residue, PTM and ion type masses from
  `res/base_data` (AME2020 masses) and the proton from `mass_constant`;
  calib.py has its own table (M 4.4e-7 Da lighter, from the sulphur mass).
  The fragment errors therefore differ by about 1e-4 ppm, and so do the
  MS/MS coefficients (480: c0 0.47366 against 0.47349 ppm); with
  TopPIC's masses put into calib.py the two fits agree to 1e-14. The
  proteoform string is parsed for single-residue fixed modifications
  (`(C)[Carbamidomethyl]`, found with `PtmBase::containsAbbrName`) and an
  `[Acetyl]-` prefix; any other PrSM gives no fragment calibrants.
- **Results** (toppic_test experiment 152, amendment 2; 9b4adfda):
  `topcalib apply` with calib.py's models writes recal_msalign.py's files
  byte for byte (480_F3_02 4,323 spectra, Bruker 4,337); `topcalib fit` on
  the zero-shift searches has calib.py's calibrant numbers (480 1,252
  precursors, 51,894 fragments of 1,381 PrSMs; Bruker 224 / 3,657 / 238)
  and MS1 coefficients to 1e-14; the corrected files of the C++ chain
  differ from the Python chain's by at most 5e-6 Da (0.006 ppm) in 40% /
  53% of the fragment lines, precursor lines identical. toppic main
  915bb387 at `-e 5` on the C++-corrected files: 480_F3_02 setting d 3,493
  PrSMs / 2,090 proteoforms / 448 proteins against the Python chain's
  3,492 / 2,089 / 448, Bruker 848 / 115 / 79 against 851 / 117 / 79 (the
  same targets at E <= 1e-5; the differences are spectra near the cutoff).
  fit 0.24 s (calib.py 1.7 s), apply 0.7 s (recal_msalign.py 1.1 s).
- **Built into toppic: the hidden `--mass-calib`** (argument `massCalib`;
  toppic_test experiment 152, amendment 3). `TopPICProgress_multi_file`
  runs `TopPIC_mass_calib` for every file, after the `-g` conversion and
  before the searches:
    - **The calibration search** is `TopPIC_identify` on a copy of the
      arguments: `shiftNumber` 0, no variable PTM file, no tag filter,
      TopMG disabled, no deamidation selector, and precursor and fragment
      tolerances of at least `min_tole` (10 ppm), so that an offset of a few
      ppm does not cost the calibrants their fragments and the index of the
      default tolerance serves (the index file names carry the fragment
      tolerance as the arguments write it, "10"). Its parameter block is
      printed as for any search; its files are deleted by `cleanToppicDir`
      with `copy_outputs` false (kept with `-k`).
    - **The calibrants** (`prsm/prsm_mass_calib::readCalibrants`, on
      `toppic_pre_raw_prsm` read with its spectra): target PrSMs (no
      `DECOY_`) with E-value <= 1e-5, no unexpected shift and no variable
      PTM; the precursor error of `getOriPrecMass()` against the proteoform
      mass, the fragment errors of the deconvoluted masses against
      `theo_peak_util::geneProteoformTheoPeak` (the activation's ion types,
      `SpPara::getMinMass`, the N-terminal label), both through
      `mass_calib_fit::addPrecCalibrant` / `addFragCalibrants`, which
      topcalib's table reader uses too (the `ms` layer cannot include
      `prsm`, hence the split). With fewer than 3 precursor or 4 fragment
      calibrants a warning is printed and the file is searched uncorrected.
    - **The corrected spectra**: the model `<base>_mass_calib.json`, the
      spectra `<base>_calib_ms2.msalign` (`msalign_mass_calib::process`) and
      a copy of the feature file; the file list then names the corrected
      file (with `-C` the earlier run's, when it exists). Both stay after the
      run, as `-g`'s `_mono` files do.
    - **The post mass matching** (`TopPIC_postprocess`, from the per-file
      argument `massCalibModelFile`) passes the model to
      `prsm_post_mass_match::process`, which corrects every raw peak's m/z
      with the MS/MS level at the spectrum's retention time (`MsHeader`
      keeps seconds, the msalign file minutes); the database's peaks are not
      changed, and the envelopes it adds to `ms2_env` are at the corrected
      m/z. `--second-precursor`'s residual pass inherits the model through
      its copy of the arguments.
    - **Names**: `sqlFileName` and `postBaseName` strip `_calib` (after
      `_resid`, before `_mono`), so with post matching the result files are
      `<base>_post_ms2_*` as without the option (a run without it in the
      same directory is replaced), with `-D` `<base>_calib_ms2_*`; `-c`
      strips `_calib` from the fractions' record names. The parameter block
      prints "Mass calibration: two-pass (zero-shift calibration search)"
      only with the option.
    - **Research switches** in `TOPPIC_MASS_CALIB_PARA` ("name=value,..."):
      `max_evalue` (1e-5), `tol_ppm` (20), `min_tole` (10), `post_match`
      (1; 0: the post mass matching reads the peaks uncorrected).
    - Without the option the run writes main 2863fa08's files (a
      500-spectrum subset of 480_F3_02, setting d: every file but the
      tables' time lines, every SQLite table).
    - **Results** (toppic_test experiment 152, amendment 3; e69a52c4,
      `-u 16`): 480_F3_02 setting d with `--mass-calib -e 5
      --disable-rescore` 3,492 PrSMs / 2,089 proteoforms / 448 proteins,
      as the Python chain of experiment 152 (topcalib's corrected file
      3,493 / 2,090 / 448, no option at 10 ppm 3,374 / 2,004 / 440); with the
      default rescoring 3,667 / 2,103 / 462 against 3,573 / 2,041 / 447.
      Bruker (report 0151's line-filtered TopFD output, `-d -t FDR -T FDR -y
      FDR`, `--disable-rescore`) 852 / 116 / 79 against topcalib's 848 / 115
      / 79 and 806 / 83 / 72 without the option. The built-in models equal
      topcalib's to 2e-3 ppm in c0 (the same 1,252 / 224 precursor
      calibrants; 51,905 / 3,661 fragment calibrants against 51,894 /
      3,657, the difference not examined; the Bruker MS/MS slope 0.315
      against 0.333 ppm per 1,000 m/z). The post mass matching on corrected
      peaks adds 11,533 masses on Bruker (offset -4 ppm) against 10,186 on
      uncorrected ones, and 503 reported PrSMs gain matched fragments (203
      lose, +934 in all; the identifications are equal); on 480 (+0.5 ppm)
      22,388 against 22,412, no change. The calibration search takes 1.6
      min (480) and 2.6 min (Bruker); with `-e 5` the 480 run takes as long
      as without the option at 10 ppm (11.3 against 11.2 min).

## Source layout

`src` is the include root, so headers are included by their path from `src`
(`common/...`, `seq/...`, `ms/spec/...`, `para/...`, `sql/...`). All of the
following compile into the single `toppic_common` **OBJECT** library (the
`COMMON_SRCS` glob in `CMakeLists.txt` lists each directory) — its object files
are linked directly into each executable, so **no `libtoppic_common.so`/`.a` is
built**; the executables are self-contained. Each layer depends only on the ones
above it, never the reverse:

- `src/common` — foundation: `base`, `util`, `xml`, `thread`, `emass`
  (isotope distributions; see its section).
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
  `factory`, `env` (envelope detection), `feature`, `mzml`, `calib` (the
  two-pass m/z recalibration of topcalib, see its section).
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
  one, and the peaks of all spectra of a longer MGF file would be merged)
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
  `prsm_second_prec` is the second precursor search (see its section);
  `prsm_mass_calib` reads the calibrants of `--mass-calib` (see the
  two-pass m/z recalibration section).
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
- `src/stat` — E-value/p-value estimation (`count`/`local`/`mcmc`/`tdgf`),
  the PTM localization (`local`, and its alternative `local_align`, see
  its section) and the rescoring of the PrSMs before the FDR (`rescore`, default with `-d`,
  see its section);
  `tdgf/cluster_prot_corrector` (see the post mass matching section) is
  here rather than in `prsm` because it recomputes E-values.
- `src/merge` — the topdiff feature-merge backend (`feature_prsm`,
  `feature_sample_merge`), used by `console/topdiff_process`.
- `src/toplib` — toplib_build's spectral library building and
  toplib_search's library search (see their sections); it depends only on
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
