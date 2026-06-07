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

project(TopPIC_claude VERSION 1.0.0 LANGUAGES C CXX)
message(STATUS "Project name: " ${PROJECT_NAME})
message(STATUS "System name: " ${CMAKE_SYSTEM_NAME})

# C++ standard (replaces a raw -std=c++17 flag).
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to an optimized build when the user did not pick a configuration.
# (Previously -O3 was hard-coded into CMAKE_CXX_FLAGS, which fought the build type.)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../bin)
message(STATUS "Executable output directory: " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# C++ flags
# -Wall is used to enable all compiler warnings
# -Wno-unknown-pragmas is used to disable pragmas warnings
# -DWITHOUT_MZ5 is used in pwiz to exclude code for mz5 format support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -DWITHOUT_MZ5")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTOPPIC_LINUX")
  # --disable-new-dtags is needed to first search the shared library in the directory ../lib/toppic
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath=$ORIGIN/../lib/toppic,--disable-new-dtags")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # fpermissive is needed to compile ext/pwiz/utility/misc/Filesystem.cpp
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -DTOPPIC_WINDOWS -DWIN32")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

message(STATUS "Cmake CXX flags: " ${CMAKE_CXX_FLAGS})
message(STATUS "Cmake shared linker flags: " ${CMAKE_SHARED_LINKER_FLAGS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# code shared by several executables
# CONFIGURE_DEPENDS makes CMake re-glob (and so reconfigure) when source files
# are added or removed under src/common.
file(GLOB_RECURSE COMMON_SRCS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp)

add_library(toppic_common SHARED ${COMMON_SRCS})
target_include_directories(toppic_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Install the shared library using the prefix-relative GNUInstallDirs locations.
# On Linux/macOS it is a LIBRARY artifact and lands in <prefix>/lib/toppic
# (matched by the $ORIGIN/../lib/toppic rpath); on Windows the DLL is a RUNTIME
# artifact (placed in <prefix>/bin next to the executables) and the import
# library is an ARCHIVE artifact.
install(TARGETS toppic_common
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
