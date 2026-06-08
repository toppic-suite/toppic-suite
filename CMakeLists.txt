cmake_minimum_required(VERSION 3.16)

# Default to clang / clang++ unless the user already chose a compiler.
# This must happen before project(), because the compiler cannot be changed
# afterwards. We only set it when unset so the user can still override with
# -DCMAKE_CXX_COMPILER=... or a toolchain file.
if(NOT CMAKE_C_COMPILER)
  find_program(CLANG_EXECUTABLE clang)
  if(CLANG_EXECUTABLE)
    set(CMAKE_C_COMPILER "${CLANG_EXECUTABLE}" CACHE STRING "C compiler")
  endif()
endif()
if(NOT CMAKE_CXX_COMPILER)
  find_program(CLANGXX_EXECUTABLE clang++)
  if(CLANGXX_EXECUTABLE)
    set(CMAKE_CXX_COMPILER "${CLANGXX_EXECUTABLE}" CACHE STRING "C++ compiler")
  endif()
endif()

project(TopPIC_claude VERSION 1.9.0 LANGUAGES C CXX)
message(STATUS "Project name: " ${PROJECT_NAME})
message(STATUS "System name: " ${CMAKE_SYSTEM_NAME})

# Version is defined once, here, via project(). src/common/util/version.hpp is
# generated from version.hpp.in so the runtime Version::getVersion() string and
# CMake's PROJECT_VERSION (used for packaging) can never drift apart.
set(TOPPIC_VERSION_SUFFIX "_alpha")   # pre-release tag appended to PROJECT_VERSION
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/util/version.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/gen/common/util/version.hpp
  @ONLY)
message(STATUS "Project version: " ${PROJECT_VERSION}${TOPPIC_VERSION_SUFFIX})

# C++ standard (replaces a raw -std=c++17 flag).
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to an optimized build when the user did not pick a configuration.
# (Previously -O3 was hard-coded into CMAKE_CXX_FLAGS, which fought the build type.)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Export compile_commands.json (used by include-what-you-use's iwyu_tool.py,
# clangd, and clang-tidy).
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optional include-what-you-use (IWYU) analysis. Off by default because IWYU is
# slow, noisy, and not always installed. Enable with -DTOPPIC_ENABLE_IWYU=ON
# (requires the include-what-you-use tool, e.g. `apt install iwyu`). IWYU runs
# alongside the compiler and prints its add/remove-include suggestions as build
# output; the build itself still proceeds normally.
option(TOPPIC_ENABLE_IWYU "Run include-what-you-use during compilation" OFF)
if(TOPPIC_ENABLE_IWYU)
  find_program(IWYU_EXECUTABLE NAMES include-what-you-use iwyu)
  if(IWYU_EXECUTABLE)
    message(STATUS "include-what-you-use enabled: ${IWYU_EXECUTABLE}")
    # Pass options through to IWYU itself with -Xiwyu (mapping files, verbosity,
    # etc.); none are needed for now.
    set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE "${IWYU_EXECUTABLE}")
  else()
    message(WARNING
      "TOPPIC_ENABLE_IWYU=ON but include-what-you-use was not found; "
      "install it (e.g. 'apt install iwyu') or pass -DIWYU_EXECUTABLE=/path/to/it.")
  endif()
endif()

# Homebrew on Apple Silicon installs to /opt/homebrew (Intel Macs use /usr/local).
# Neither prefix is on the default compiler / linker search path, so the project's
# bare-name target_link_libraries calls and #include directives need help here.
# NOTE: this must come AFTER project(), because the APPLE variable is only set by project().
if(APPLE)
  include_directories(/opt/homebrew/include)
  link_directories(/opt/homebrew/lib)
endif()

# Install directories, expressed RELATIVE to CMAKE_INSTALL_PREFIX (default
# /usr/local on Unix) so that a custom prefix, DESTDIR staging and packaging all
# work. The previous absolute paths (e.g. "/usr/local/lib/toppic") ignored the
# prefix entirely and could not be staged or relocated.
include(GNUInstallDirs)   # standard bin/lib/share defaults, all under the prefix

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # The shared library goes in <prefix>/lib/toppic and the data in
  # <prefix>/share/toppic. NOTE: "lib/toppic" must stay in sync with the
  # "$ORIGIN/../lib/toppic" rpath set in the Linux linker flags below.
  set(CMAKE_INSTALL_LIBDIR  "lib/toppic")
  set(CMAKE_INSTALL_DATADIR "${CMAKE_INSTALL_DATADIR}/toppic")
ELSEIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # macOS: resources are installed next to the executable so that
  # getResourceDir()'s "<exec_dir>/resources" lookup finds them.
  set(CMAKE_INSTALL_DATADIR "${CMAKE_INSTALL_BINDIR}/resources")
ENDIF()

message(STATUS "Install prefix:              ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Bin installation directory:  ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
message(STATUS "Lib installation directory:  ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
message(STATUS "Data installation directory: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}")

# Compile the resolved (absolute) data directory into the binaries so that
# file_util::getSharedDirName() / getResourceDir() look exactly where the
# resources are installed, instead of a hard-coded /usr/share/toppic.
if(IS_ABSOLUTE "${CMAKE_INSTALL_DATADIR}")
  set(TOPPIC_SHARED_DIR "${CMAKE_INSTALL_DATADIR}")
else()
  set(TOPPIC_SHARED_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}")
endif()
add_compile_definitions(TOPPIC_SHARED_DIR="${TOPPIC_SHARED_DIR}")
message(STATUS "Compiled-in shared dir:      ${TOPPIC_SHARED_DIR}")

# Install the runtime resources (base_data, envcnn_models, ecscore_models,
# p_value_table, topmsv) into the compiled-in shared dir. The trailing slash on
# "resources/" installs the directory's CONTENTS into CMAKE_INSTALL_DATADIR, so
# at runtime getResourceDir()/getSharedDirName() (== TOPPIC_SHARED_DIR ==
# <prefix>/${CMAKE_INSTALL_DATADIR}) finds base_data/, envcnn_models/, ...
# directly inside it. On macOS CMAKE_INSTALL_DATADIR is "<bindir>/resources",
# which getResourceDir()'s "<exec_dir>/resources" lookup also resolves.
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/resources/
        DESTINATION ${CMAKE_INSTALL_DATADIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../bin)
message(STATUS "Executable output directory: " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# C++ flags
# -Wall is used to enable all compiler warnings
# -Wno-unknown-pragmas is used to disable pragmas warnings
# -DWITHOUT_MZ5 is used in pwiz to exclude code for mz5 format support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -DWITHOUT_MZ5")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTOPPIC_LINUX")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # fpermissive is needed to compile ext/pwiz/utility/misc/Filesystem.cpp
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -DTOPPIC_WINDOWS -DWIN32")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

message(STATUS "Cmake CXX flags: " ${CMAKE_CXX_FLAGS})
message(STATUS "Cmake shared linker flags: " ${CMAKE_SHARED_LINKER_FLAGS})

# pugixml is used by src/common/xml for XML parsing and serialization.
find_package(pugixml REQUIRED)
message(STATUS "pugixml version: ${pugixml_VERSION}")

# A trimmed copy of htslib (faidx/bgzf/hfile only) is vendored under ext/htslib;
# seq/fasta_index_reader uses faidx to random-access indexed FASTA files.
# It is built as a small static library and linked into toppic_common. The
# vendored C sources include themselves as "htslib/<name>.h", so the parent
# ext/ directory is the include root. zlib (and threads, for BGZF_MT) are the
# only external dependencies.
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# SQLite3 backs sql/sql_util and ms/mzml/mzml_ms_sql_writer (the topfd SQLite
# output).
find_package(SQLite3 REQUIRED)
file(GLOB HTSLIB_SRCS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib/*.c)
add_library(htslib STATIC ${HTSLIB_SRCS})
target_include_directories(htslib SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ext)
target_link_libraries(htslib PUBLIC ZLIB::ZLIB Threads::Threads)
# The vendored C sources are third-party; keep their (many) warnings out of our
# build output and build them with position-independent code so they can be
# linked into the toppic_common shared library.
set_target_properties(htslib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(htslib PRIVATE -w)

# A trimmed copy of ProteoWizard (pwiz) is vendored under ext/pwiz; ms/mzml's
# pw_ms_reader uses it to read mzML/mzMLb files. Only the four source
# directories below are compiled (the rest of the tree is headers). It is built
# as a static library against Boost + zlib and linked into toppic_common. The
# vendored sources include themselves as "pwiz/...", so ext/ is the include root.
# WITHOUT_MZ5 drops the HDF5-backed mz5 format so no HDF5 dependency is needed.
find_package(Boost 1.74 REQUIRED COMPONENTS filesystem iostreams thread chrono system serialization program_options)
file(GLOB PWIZ_SRCS CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/utility/minimxml/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/utility/misc/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/data/common/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/data/msdata/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/data/msdata/mzmlb/*.cpp)
add_library(pwiz STATIC ${PWIZ_SRCS})
target_include_directories(pwiz SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ext)
target_link_libraries(pwiz PUBLIC
  Boost::filesystem Boost::iostreams Boost::thread Boost::chrono Boost::system
  ZLIB::ZLIB)
target_compile_definitions(pwiz PUBLIC WITHOUT_MZ5)
set_target_properties(pwiz PROPERTIES POSITION_INDEPENDENT_CODE ON)
# Third-party: silence its warnings (the -Wdeprecated noise from Boost especially).
target_compile_options(pwiz PRIVATE -w)

# src for the hand-written headers, and the build-tree gen/ dir for the
# CMake-generated common/util/version.hpp.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/gen)

# Code shared by every executable. It is compiled once into an OBJECT library
# (no libtoppic_common.so / .a is produced); the object files are linked
# directly into each command-line / GUI executable, which keeps them
# self-contained (no shared library to locate at runtime). Its PUBLIC/PRIVATE
# link dependencies and usage requirements still propagate to the executables
# that link it.
# CONFIGURE_DEPENDS makes CMake re-glob (and so reconfigure) when source files
# are added or removed under these directories.
file(GLOB_RECURSE COMMON_SRCS CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/seq/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/para/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/sql/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/topfd/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/topdia/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/prsm/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/search/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/stat/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/merge/*.cpp)

add_library(toppic_common OBJECT ${COMMON_SRCS})
target_include_directories(toppic_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_BINARY_DIR}/gen)

# PUBLIC: the common/xml headers expose pugixml types, so consumers of
# toppic_common need pugixml on their include and link lines too.
target_link_libraries(toppic_common PUBLIC pugixml::pugixml)

# seq/fasta_index_reader.hpp includes <htslib/faidx.h>, so the vendored
# htslib include path is PUBLIC; the library itself is an implementation detail
# (PRIVATE) of toppic_common.
target_include_directories(toppic_common SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ext)
target_link_libraries(toppic_common PRIVATE htslib)

# sql/sql_util.hpp and ms/mzml/mzml_ms_sql_writer.hpp expose <sqlite3.h>, so the
# SQLite3 include path is PUBLIC; the library is PRIVATE (an implementation
# detail of toppic_common).
target_link_libraries(toppic_common PUBLIC SQLite::SQLite3)

# ms/mzml/pw_ms_reader.hpp exposes pwiz types (pwiz::msdata::*), so consumers of
# toppic_common need the pwiz headers (and Boost) too: link pwiz PUBLIC.
target_link_libraries(toppic_common PUBLIC pwiz)

# filter/massmatch/mass_match (de)serializes its index with Boost.Serialization;
# its header exposes boost::serialization::access, so link it PUBLIC.
target_link_libraries(toppic_common PUBLIC Boost::serialization)

# ONNX Runtime: topfd/envcnn and topfd/ecscore/score run the EnvCNN / ECScore
# neural-network models via the ONNX Runtime C++ API. A prebuilt libonnxruntime
# and the C/C++ API headers are vendored under ext/onnx (the sources include it
# as "onnx/onnxruntime_cxx_api.h", so the existing SYSTEM ext/ include covers the
# headers). The library is an IMPORTED shared object linked into toppic_common.
add_library(onnxruntime SHARED IMPORTED)
set_target_properties(onnxruntime PROPERTIES
  IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/ext/onnx/libonnxruntime.so.1.14.1
  IMPORTED_SONAME libonnxruntime.so.1.14.1)
target_link_libraries(toppic_common PRIVATE onnxruntime)

# ---------------------------------------------------------------------------
# Command-line executables
#
# Each tool is its own main() in src/console plus its argument parser and, for
# the four multi-step tools, its src/console/<tool>_process.cpp orchestrator
# (topfd/topdia drive the topfd_process/topdia_process orchestrators that are
# already in toppic_common). Everything else they call lives in the library, so
# each exe links toppic_common + Boost::program_options (the argument parsers).
# ---------------------------------------------------------------------------
function(toppic_console_exe name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE toppic_common Boost::program_options)
endfunction()

toppic_console_exe(topfd
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_argument.cpp)
toppic_console_exe(topdia
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdia.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdia_argument.cpp)
toppic_console_exe(topindex
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topindex.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topindex_argument.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topindex_process.cpp)
toppic_console_exe(toppic
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_argument.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_process.cpp)
toppic_console_exe(topmg
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_argument.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_process.cpp)
toppic_console_exe(topdiff
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff_argument.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff_process.cpp)

# ---------------------------------------------------------------------------
# GUI executables (Qt5)
#
# Each tool's dialog collects parameters and shells out to the matching
# command-line tool via QProcess (gui/util/command). All the heavy lifting
# (common/seq/para/ms/prsm/filter/search/stat/topfd/topdia) already
# lives in the toppic_common shared library, so a GUI target only compiles its
# own Qt sources -- plus, for the five tools whose dialog reads default values
# from the console parser, that tool's src/console/<tool>_argument.cpp -- and
# links against toppic_common, Qt5 and Boost::program_options.
#
# (The src/gui/topmerge sources are migrated for completeness but, as upstream,
# get no executable target: the merge backend is not part of this library.)
# ---------------------------------------------------------------------------
find_package(Qt5 COMPONENTS Widgets Core Gui REQUIRED)

# Shared by every GUI target (a QProcess command builder + message helpers).
# It is QObject-based, so it is built per target via the target's AUTOMOC.
file(GLOB GUI_UTIL_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/util/*.cpp)

# toppic_gui_exe(<target> <gui-subdir> [extra sources...])
function(toppic_gui_exe name dir)
  file(GLOB _srcs
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/${dir}/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/${dir}/*.qrc)
  add_executable(${name} ${_srcs} ${GUI_UTIL_SRCS} ${ARGN})
  # AUTOUIC resolves the dialogs' "gui/<dir>/ui_*.h" includes by finding the
  # matching .ui under src/ (the search path) and emitting the header with the
  # same relative path into the target's autogen include dir.
  set_target_properties(${name} PROPERTIES
    AUTOMOC ON AUTOUIC ON AUTORCC ON
    AUTOUIC_SEARCH_PATHS ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(${name} PRIVATE
    toppic_common Qt5::Widgets Qt5::Core Qt5::Gui Boost::program_options)
endfunction()

toppic_gui_exe(topfd_gui    topfd)
toppic_gui_exe(topindex_gui topindex ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topindex_argument.cpp)
toppic_gui_exe(toppic_gui   toppic   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_argument.cpp)
toppic_gui_exe(topmg_gui    topmg    ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_argument.cpp)
toppic_gui_exe(topdiff_gui  topdiff  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff_argument.cpp)
toppic_gui_exe(topdia_gui   topdia   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdia_argument.cpp)
