# Vendored third-party code

This directory holds trimmed copies of the third-party libraries that TopPIC
Suite builds or links directly. They are third-party code: do not reformat
them or hold them to the project's include-order and IWYU rules.

| Directory | Contents | Version / source |
|---|---|---|
| `boost/` | The Boost headers that `pwiz/` includes and the system Boost does not provide in the version pwiz expects: the `boost::enums` library, the singleton utility and Boost.Nowide (`nowide/system.hpp` and `nowide/iostream.cpp` from ProteoWizard's `libraries/boost_aux`, the other `nowide` headers from the standalone Boost.Nowide that earlier ProteoWizard releases bundled). | Bundled with ProteoWizard |
| `htslib/` | Library for creating FASTA indexes and reading FASTA files (`faidx`, `bgzf`, `hfile` only). | 2016 release, <http://www.htslib.org/download/> |
| `pwiz/` | The ProteoWizard code that TopPIC uses to read mzML and mzXML files: only the files TopPIC needs (see below). | 3.0.26278 (git 0c7a916, 2026-10-05), source archive `pwiz-src-without-v-3_0_26278_0c7a916.tar.bz2`, <http://proteowizard.sourceforge.net/download.html> |
| `onnx/` | Microsoft ONNX Runtime library (the C/C++ API headers and the Linux x86-64 shared library `libonnxruntime.so.1.23.2` with its SONAME symlink `libonnxruntime.so.1`, from the `onnxruntime-linux-x64-1.23.2.tgz` release archive). | 1.23.2, <https://github.com/microsoft/onnxruntime/releases/tag/v1.23.2> |

## ProteoWizard (`pwiz/` and `boost/`)

TopPIC includes pwiz only in `src/ms/mzml/pw_ms_reader.{hpp,cpp}`. `pwiz/` and
`boost/` hold just the files that this reader needs (105 and 19 files):

- the 39 source files of pwiz that the linker takes into the TopPIC
  executables, plus `boost/nowide/iostream.cpp` (the UTF-8 console streams
  that pwiz writes its messages to on Windows; an empty translation unit
  elsewhere);
- every header that these sources and the reader reach through `#include`,
  on every platform, apart from the code that `WITHOUT_MZ5` and
  `WITHOUT_MZMLB` switch off.

`CMakeLists.txt` compiles every `.cpp` file under `pwiz/`, and
`boost/nowide/iostream.cpp`, into the static library `pwiz` with
`WITHOUT_MZ5` and `WITHOUT_MZMLB` (so no HDF5 is needed: the mz5 and mzMLb
formats are not supported, and pwiz reports an mzMLb file as such), and links
it against Boost, zlib and zstd (since 3.0.26278 pwiz's binary data encoder
decodes zstd-compressed data arrays, `pwiz/utility/misc/mzd.hpp`).

`pwiz/data/msdata/Version.cpp` is not in the source archive: pwiz's build
generates it (`libraries/generate-version.jam`). The copy here has the text
that the build generates for version 3.0.26278, revision 0c7a916, branch
master.

Local changes to the upstream files (keep them when updating):

- `pwiz/utility/misc/Filesystem.cpp`, for the MinGW (GCC) build on Windows:
  the Windows code is selected by `WIN32` instead of `_MSC_VER`; the struct
  `SYSTEM_HANDLE_INFORMATION` is renamed `PWIZ_SYSTEM_HANDLE_INFORMATION`
  (to stay clear of the MinGW Windows headers); and `GetLibraryProcAddress`
  returns `nullptr` instead of the `GetProcAddress` result, so pwiz finds no
  `ntdll.dll` functions: closing the handles of a locked file reports an
  error and is skipped, and Wine is never detected.
- `boost/nowide/iostream.cpp`: without `#include <boost/locale/utf.hpp>`.
- `boost/nowide/system.hpp`: `wcmd.c_str()` instead of `wcmd.get()`, the
  interface of the older `nowide/stackstring.hpp` vendored here.

### Updating pwiz

1. Download the source archive of the new version
   (`pwiz-src-without-v-<version>.tar.bz2`) and extract it.
2. Copy the new version of every file in `pwiz/` (from `pwiz/` of the
   archive) and `boost/` (from `libraries/boost_aux/boost/`;
   `nowide/iostream.cpp` from `libraries/boost_aux/libs/nowide/src/`; the
   older `nowide` headers that `boost_aux` no longer has stay as they are).
   Re-apply the local changes and write `pwiz/data/msdata/Version.cpp` for the
   new version.
3. Build. A missing header means that the new code includes a file that is
   not here yet: copy it, and the headers it includes. An undefined
   reference means that a symbol moved to a source file that is not here:
   copy that `.cpp` file. Files that nothing includes any more can be
   removed.

The file set of 3.0.26278 was computed this way from the archive: the
include closure of the reader's pwiz headers, a build with the `.cpp` file
of every header reached, and the linker map of `topfd`, which lists the
library members the executables use (toppic_test experiment 134).
