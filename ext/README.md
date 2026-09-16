# Vendored third-party code

This directory holds trimmed copies of the third-party libraries that TopPIC
Suite builds or links directly. They are third-party code: do not reformat
them or hold them to the project's include-order and IWYU rules.

| Directory | Contents | Version / source |
|---|---|---|
| `boost/` | Some code from Boost used by `pwiz/libraries/boost_aux/` (`nowide`, the `boost::enums` library, `foreach_field.hpp`). | Bundled with ProteoWizard |
| `htslib/` | Library for creating FASTA indexes and reading FASTA files (`faidx`, `bgzf`, `hfile` only). | 2016 release, <http://www.htslib.org/download/> |
| `pwiz/` | Some code from ProteoWizard for reading mzML and mzXML files. | 3.0_25117 (2025-04-29), <http://proteowizard.sourceforge.net/download.html> |
| `onnx/` | Microsoft ONNX Runtime library (C/C++ API headers and the Linux x86-64 shared library). | 1.14.1, <https://github.com/microsoft/onnxruntime/releases/tag/v1.14.1> |
