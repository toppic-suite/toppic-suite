#Specify the version being used as well as the language
cmake_minimum_required(VERSION 3.1)
#Name your project here
project(TopPIC)

SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../lib)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O3 -DWITHOUT_MZ5 -fPIC")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O2 -pg -fPIC")
message(STATUS "system name " ${CMAKE_SYSTEM_NAME})

# for Qt
set(CMAKE_AUTOMOC ON)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
find_package(Qt5Widgets)
find_package(Qt5Core)
find_package(Qt5Gui)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Linux specific code
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-rpath,./toppic_resources/lib")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
                      ${Qt5Widgets_INCLUDE_DIRS}
                      ${Qt5Core_INCLUDE_DIRS}
                      ${Qt5Gui_INCLUDE_DIRS})
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # Windows specific code
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWIN32")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
                      ${Qt5Widgets_INCLUDE_DIRS}
                      ${Qt5Core_INCLUDE_DIRS}
                      ${Qt5Gui_INCLUDE_DIRS})
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

message(STATUS "cmake cxx flags " ${CMAKE_CXX_FLAGS})

# pwiz
file(GLOB PWIZ_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/pwiz/utility/minimxml/*.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/src/pwiz/utility/misc/*.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/src/pwiz/data/common/*.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/src/pwiz/data/msdata/*.cpp)

# common
file(GLOB HTS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/htslib/*.c)
file(GLOB BASE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/base/*.cpp)
file(GLOB SPEC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/spec/*.cpp)
file(GLOB PRSM_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/prsm/*.cpp)
file(GLOB ZERO_PTM_FILTER  ${CMAKE_CURRENT_SOURCE_DIR}/src/zeroptmfilter/*.cpp)
file(GLOB ZERO_PTM_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/zeroptmsearch/*.cpp)
file(GLOB DIAG_FILTER_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/diagfilter/*.cpp)
file(GLOB ONE_PTM_SEARCH_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/oneptmsearch/*.cpp)
file(GLOB ONE_PTM_FILTER_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/oneptmfilter/*.cpp)
file(GLOB PTM_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/ptmsearch/*.cpp)
file(GLOB TDGF_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/tdgf/*.cpp)
file(GLOB PRSM_VIEW_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/prsmview/*.cpp)
file(GLOB LOCAL_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/local/*.cpp)
file(GLOB SUFFIX_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/suffix/*.cpp)

# toppic
file(GLOB TOPPIC_ARGU_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_argument.cpp)
file(GLOB TOPPIC_PROCESS_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_process.cpp)
file(GLOB TOPPIC_MAIN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic.cpp)
file(GLOB TOPPIC_MATCH_PEAK_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_match_peak.cpp)

# mass graph
file(GLOB GRAPH_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/graph/*.cpp)
file(GLOB GRAPH_ALIGN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/graphalign/*.cpp)

# topmg
file(GLOB TOPMG_ARGU_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_argument.cpp)
file(GLOB TOPMG_PROCESS_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_process.cpp)
file(GLOB TOPMG_MAIN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg.cpp)

# feature detection
file(GLOB FEATURE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/feature/*.cpp)

# topfd
file(GLOB TOPFD_ARGU_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_argument.cpp)
file(GLOB TOPFD_MAIN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd.cpp)

file(GLOB TOPFD_TRAIN_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_train.cpp)

# GUI
file(GLOB TOPFD_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topfd/*.cpp)
file(GLOB TOPFD_GUI_QRC   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topfd/*.qrc)
file(GLOB TOPPIC_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/toppic/*.cpp)
file(GLOB TOPPIC_GUI_QRC   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/toppic/*.qrc)
file(GLOB TOPMG_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmg/*.cpp)
file(GLOB TOPMG_GUI_QRC   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmg/*.qrc)

# MCMC
file(GLOB MCMC_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/mcmc/*.cpp)

# unit test using catch
file(GLOB UNIT_TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/*.cpp)

# quant
file(GLOB TOPQUANT_MAIN_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topquant.cpp)
file(GLOB TOPQUANT_ARGU_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topquant_argument.cpp)

add_executable(toppic ${TOPPIC_MAIN_SRCS} ${TOPPIC_PROCESS_SRCS} ${TOPPIC_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
    ${ZERO_PTM_SRCS} ${ZERO_PTM_FILTER} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${ONE_PTM_SEARCH_SRCS}
    ${PTM_SRCS} ${TDGF_SRCS} ${LOCAL_SRCS} ${PRSM_VIEW_SRCS})

add_executable(toppic_gui ${TOPPIC_GUI_SRCS} ${TOPPIC_GUI_QRC} ${TOPPIC_PROCESS_SRCS} ${TOPPIC_ARGU_SRCS} 
    ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
    ${ZERO_PTM_SRCS} ${ZERO_PTM_FILTER} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${ONE_PTM_SEARCH_SRCS}
    ${PTM_SRCS} ${TDGF_SRCS} ${LOCAL_SRCS} ${PRSM_VIEW_SRCS})

add_executable(topmg ${TOPMG_MAIN_SRC} ${TOPMG_PROCESS_SRCS} ${TOPMG_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
    ${ZERO_PTM_FILTER} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${ZERO_PTM_SRCS} ${GRAPH_SRCS}
    ${GRAPH_ALIGN_SRCS} ${TDGF_SRCS} ${MCMC_SRCS} ${ONE_PTM_SEARCH_SRCS} ${PRSM_VIEW_SRCS})

add_executable(topmg_gui ${TOPMG_GUI_SRCS} ${TOPMG_GUI_QRC} ${TOPMG_PROCESS_SRCS} ${TOPMG_ARGU_SRCS} 
    ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
    ${ZERO_PTM_FILTER} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${ZERO_PTM_SRCS} ${GRAPH_SRCS}
    ${GRAPH_ALIGN_SRCS} ${TDGF_SRCS} ${MCMC_SRCS} ${ONE_PTM_SEARCH_SRCS} ${PRSM_VIEW_SRCS})

add_executable(topfd ${TOPFD_MAIN_SRCS} ${TOPFD_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS} ${PWIZ_SRC})

add_executable(topfd_gui ${TOPFD_GUI_SRCS} ${TOPFD_GUI_QRC} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS} ${PWIZ_SRC})

add_executable(topquant ${TOPQUANT_MAIN_SRCS} ${TOPQUANT_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries(toppic xalan-c xalanMsg xerces-c
        boost_filesystem-mt boost_system-mt 
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

    target_link_libraries(topmg xalan-c xalanMsg xerces-c
        boost_filesystem-mt boost_system-mt
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

    target_link_libraries(topfd xalan-c xalanMsg xerces-c
        boost_program_options-mt boost_filesystem-mt
        boost_system-mt boost_thread-mt
        boost_iostreams-mt boost_chrono-mt pthread z Ws2_32)

    # for Qt
    SET(QTOPTION -Wl,-subsystem,windows)
    target_link_libraries(toppic_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui xalan-c xalanMsg xerces-c
        boost_filesystem-mt boost_system-mt 
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

    target_link_libraries(topfd_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui
        xalan-c xerces-c boost_filesystem-mt
		boost_system-mt boost_thread-mt 
		boost_iostreams-mt boost_chrono-mt pthread z Ws2_32)

    target_link_libraries(topmg_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui xalan-c xalanMsg xerces-c
        boost_filesystem-mt boost_system-mt
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

    target_link_libraries(topquant xalan-c xerces-c boost_program_options-mt
        boost_filesystem-mt boost_system-mt boost_thread-mt pthread z Ws2_32)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_executable(toppic-test ${UNIT_TEST_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SUFFIX_SRCS})

    add_executable(topfd_train ${TOPFD_TRAIN_SRCS} ${TOPFD_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS} ${PWIZ_SRC})

    add_executable(toppic_match_peak ${TOPPIC_MATCH_PEAK_SRCS} ${TOPPIC_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS} ${PRSM_VIEW_SRCS})

    target_link_libraries(toppic xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(toppic_gui Qt5Widgets Qt5Core Qt5Gui xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(toppic_match_peak xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topmg xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topfd xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread boost_iostreams boost_chrono pthread z)

    target_link_libraries(topfd_train xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread boost_iostreams boost_chrono pthread z)

    target_link_libraries(topfd_gui Qt5Widgets Qt5Core Qt5Gui
        xalan-c xerces-c boost_filesystem boost_system boost_thread boost_iostreams boost_chrono pthread z)

    target_link_libraries(topmg_gui Qt5Widgets Qt5Core Qt5Gui xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(toppic-test xalan-c xerces-c boost_filesystem
        boost_system pthread z)

    target_link_libraries(topquant xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

