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
  chrono system serialization program_options)`. `ms/mzml`'s pwiz reader needs
  `filesystem`/`iostreams`/`thread`/`chrono`/`system`; `filter/massmatch/
  mass_match` needs `serialization` (its binary index archive — linked PUBLIC
  into `toppic_common` because the header exposes `boost::serialization::access`);
  the GUI argument parsers need `program_options`. `search/graph` additionally
  uses the **header-only** Boost Graph Library (no link). Prefer `std::mutex` etc.
  over `boost::*` in our own code (the filter/index mutexes were converted from
  `boost::mutex`; only `mass_match`'s serialization genuinely needs Boost). (Boost
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
  guard; `docs/build/windows_build.md` tells users to install `...-gcc`, not `-clang`.
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
  `libonnxruntime.so.1.14.1` and the C/C++ API headers are vendored under
  `ext/onnx` (sources include it as `"onnx/onnxruntime_cxx_api.h"`, covered by
  the `SYSTEM` `ext/` include); it is linked as an `IMPORTED` shared object,
  `PRIVATE`, into `toppic_common`. The vendored `.so` is Linux x86-64 only;
  on macOS and Windows `CMakeLists.txt` instead locates the system package
  (Homebrew's `libonnxruntime.dylib`; MSYS2 UCRT64's `libonnxruntime.dll.a`)
  with `find_library` (overridable with `-DONNXRUNTIME_LIBRARY=...`) as an
  `UNKNOWN IMPORTED` target, and still compiles against the vendored 1.14
  headers. Windows also links `ws2_32` into the vendored htslib and pwiz
  (Winsock), builds the GUI targets as `WIN32_EXECUTABLE`, and installs `res/`
  next to the executables because `getResourceDir()` does not consult the
  compiled-in shared dir there. Treat it as third-party.

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
- On Linux the vendored `ext/onnx/libonnxruntime.so.1.14.1` (which every
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

The build instructions are split by platform. `README.md` carries the full
Ubuntu walkthrough (packages, clone with Git LFS, configure/build, install,
uninstall) and only a one-paragraph link for each other platform; the details
live in `docs/build/`, all laid out with the same six steps:

- `docs/build/windows_build.md` — MSYS2 **UCRT64** shell, Ninja generator, GCC (not
  clang; see the compiler-default note above), `pacman` package list, system
  ONNX Runtime, `res/` installed next to the executables.
- `docs/build/redhat_build.md` — RHEL 10 (and Rocky/AlmaLinux): CodeReady Builder +
  EPEL for `pugixml-devel`, `dnf` package list, otherwise identical to Ubuntu.
- `docs/build/macos_build.md` — Xcode Command Line Tools + Homebrew, Homebrew
  `onnxruntime` (nothing is vendored for macOS; `-DONNXRUNTIME_LIBRARY`
  override), Qt6 hint, resources installed to `<prefix>/bin/res`.

When you change anything platform-specific in `CMakeLists.txt` — a required
package, the compiler choice, how a library is located, the install layout —
update the matching doc (and the README for Ubuntu) in the same commit; the
docs are hand-written and nothing checks them against the build.

## TopFD user manual (`docs/topfd_manual.md`)

`docs/topfd_manual.md` is the user-facing manual for the `topfd` command line:
its two functions (mzML/mzXML deconvolution, and `-T` text-peak-list
deconvolution of a single MS/MS spectrum), the input requirements, every
output file, and every option with its default — the visible ones in the
tables of section 1.4 and the hidden ones (`-k`, `-M`, `-O`,
`--max-miss-peak-num`, `--disable-filter-by-mz`, `--output-dp-envs`) in the
"Advanced options" table. It was written from `console/topfd_argument.cpp`,
`topfd/common/topfd_process.cpp`, `topfd/common/topfd_single_process.cpp` and
the `msalign`/feature/env writers, and nothing checks it against the code, so:

- When you add, rename, hide/unhide or change the default of an option in
  `topfd_argument.cpp`, or change an output file name, format or the columns
  of a writer, update the manual in the same commit. The option table there
  mirrors the `-h` text; keep the two descriptions saying the same thing.
- Do not document an option in the manual that the parser does not accept,
  and do not put a `-T`-mode claim in section 2 without checking
  `topfd_single_process.cpp` — that mode ignores the MS1/feature options.
- The examples use no `-a`: the activation defaults to `FILE` (read from the
  input), which is the recommended usage; keep new examples consistent.

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

## Post mass matching (`prsm/prsm_post_mass_match`, toppic only)

After the search and the match recount, `console/toppic_process` runs
`prsm_post_mass_match::process`, which matches the theoretical fragment
masses of each PrSM that no deconvoluted mass matched against the
**centroided MS/MS peaks stored in TopFD's SQLite database** (MSPathFinderT's
envelope test), scores the added masses with EnvCNN, and writes
`<base>_post_ms2.msalign` plus PrSM/feature files under that name, which every
downstream step (clustering, FDR, tables, SQL output) then reads instead of
the input. Facts to keep straight:

- **It is on by default.** `initArguments()` sets `postMassMatch = "true"`;
  `-E` / `--disable-post-match` sets it to false (there is no enabling flag
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
- **Result file names carry `_post_ms2`** whenever post matching ran; the
  cleanup in `toppic_process` (`cleanToppicDir`) takes that name to find the
  intermediate files. With `-E` the names are the plain `<base>_ms2` ones.
- **GUI wiring.** `toppicwindow.ui` has `postMassMatchCheckBox` (checked by
  default from the parser's argument map, enables `postMinPeakNumEdit`) on the
  bottom row of the Advanced Parameters tab; `gui/util/command.cpp` emits
  `-E` when the box is **cleared** and `-I <n>` only while it is checked.
- **Short options in `toppic_argument.cpp` are nearly exhausted.** Before
  adding one, list the letters in use with
  `grep -oE '"[a-z-]+,[A-Za-z]"' src/console/toppic_argument.cpp`; a clash is
  only reported at run time as "option '-X' is ambiguous" (that is how `-P`
  collided with `--proteoform-ppm-error`). Every option is declared twice
  (the visible `display_desc` and the full parse `desc`); change both.

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
  cutoff (`SpPara::getEnvCnnCutoff()`, default 0.2): masses whose fourth
  msalign column is below it are skipped, keeping the file's peak ids. The
  filter, search and E-value/p-value stages (`filter/*`, `search/*`,
  `stat/tdgf`, `stat/mcmc`, `search/graph/spec_graph_reader`) construct their
  readers with it, so they work on the filtered masses; the PrSM output
  stages (`prsm_reader_util`, coverage, localization, tables) deliberately
  read the unfiltered spectra. Use the cutoff constructor in a new
  identification stage, and the three-argument one in a new output stage.
  Because the search stores matched mass/fragment numbers counted on the
  filtered masses, toppic and topmg run `prsm/prsm_match_num_recount` on
  `<tool>_raw_prsm` (writing `<tool>_recount`) before clustering to recount
  them against the full spectra; the number of masses in the tables is not
  stored in the PrSM but computed from the attached full spectra at output.
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
- `src/filter` — proteoform filtering (`diag`/`index`/`massmatch`/`mng`/
  `oneptm`/`varptm`/`zeroptm`). `massmatch/mass_match` (de)serializes its binary
  index with Boost.Serialization. The file-global `serial_mutex` (and the other
  filter/index mutexes) live in **anonymous namespaces** — they previously had
  external linkage and the same name in several TUs, which only becomes a
  multiple-definition error once everything compiles into one library.
- `src/search` — sequence/spectrum alignment (`diag`/`graph`/`graphalign`/
  `oneptmsearch`/`ptmsearch`/`varptmsearch`/`zeroptmsearch`). `graph` uses the
  header-only Boost Graph Library (`adjacency_list`, `graph_traits`, `graphviz`).
- `src/stat` — E-value/p-value estimation (`count`/`local`/`mcmc`/`tdgf`).
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
