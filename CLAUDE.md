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
  `sql/sql_util` and `ms/mzml/mzml_ms_sql_writer`.
- **rapidjson** — header-only, on the default include path (no `find_package`,
  no link); used by `ms/mzml/mzml_ms_json_writer`.
- **Boost uBLAS** — header-only; used by `ms/util` (Savitzky-Golay). Its headers
  still derive from the C++17-deprecated `std::iterator`, so wrap Boost includes
  in `#pragma GCC diagnostic ignored "-Wdeprecated-declarations"` to keep the
  build warning-free. (Prefer `std::mutex` etc. over `boost::*` in our own code.)

**ProteoWizard (pwiz)** is *not* available here. `ms/mzml/pw_ms_reader` and
`ms/mzml/mzml_ms_group_reader` (which includes it) read mzML via pwiz and are
therefore **not yet migrated** — pwiz is vendored as `ext/pwiz` upstream and is
too large to bring in. Nothing else depends on them.

## Source layout

`src` is the include root, so headers are included by their path from `src`
(`common/...`, `seq/...`, `ms/spec/...`, `para/...`, `sql/...`). All of the
following compile into the single `toppic_common` shared library (the
`COMMON_SRCS` glob in `CMakeLists.txt` lists each directory). Each layer depends
only on the ones above it, never the reverse:

- `src/common` — foundation: `base`, `util`, `xml`, `thread`.
- `src/sql` — thin SQLite helper (`sql_util`).
- `src/para` — analysis parameters (`sp_para`, `peak_tolerance`).
- `src/seq` — sequence/proteoform layer.
- `src/ms` — mass-spectrum layer: `spec` (peaks/spectra/msalign), `util`,
  `msmap`, `factory`, `env` (envelope detection), `feature`, `mzml`.

When migrating a folder from the upstream Xerces tree, watch for include guards
that don't match the destination path (e.g. an `ms/env` file guarded
`TOPPIC_TOPFD_ENV_*`) and rename them to `TOPPIC_<PATH>_<FILE>_HPP_`.
