#Specify the version being used as well as the language
cmake_minimum_required(VERSION 3.1)
#Name your project here
project(TopPIC)

SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../lib)

SET(CMAKE_INSTALL_PREFIX $ENV{HOME}/toppic)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O3 -DWITHOUT_MZ5 -fPIC")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall  -g -std=c++11  -DWITHOUT_MZ5 -fPIC")
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
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTOPPIC_LINUX")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
                      ${Qt5Widgets_INCLUDE_DIRS}
                      ${Qt5Core_INCLUDE_DIRS}
                      ${Qt5Gui_INCLUDE_DIRS})
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # Windows specific code
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTOPPIC_WINDOWS -DWIN32")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
                      ${Qt5Widgets_INCLUDE_DIRS}
                      ${Qt5Core_INCLUDE_DIRS}
                      ${Qt5Gui_INCLUDE_DIRS})
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

message(STATUS "cmake cxx flags " ${CMAKE_CXX_FLAGS})

# pwiz
file(GLOB PWIZ_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/pwiz/utility/minimxml/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/pwiz/utility/misc/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/pwiz/data/common/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/pwiz/data/msdata/*.cpp)

# hts
file(GLOB HTS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/htslib/*.c)

# common
file(GLOB_RECURSE COMMON_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp)
file(GLOB SEQ_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/seq/*.cpp)
file(GLOB SPEC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/spec/*.cpp)
file(GLOB FEATURE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/feature/*.cpp)
file(GLOB PRSM_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/prsm/*.cpp)

# suffix tree
file(GLOB SUFFIX_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/suffix/*.cpp)

# test
file(GLOB TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/*.cpp)

# classes for toppic
file(GLOB_RECURSE FILTER_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/*.cpp)

file(GLOB ZERO_PTM_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/zeroptmsearch/*.cpp)
file(GLOB ONE_PTM_SEARCH_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/oneptmsearch/*.cpp)
file(GLOB PTM_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/ptmsearch/*.cpp)

file(GLOB TDGF_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/tdgf/*.cpp)
file(GLOB PRSM_VIEW_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/prsmview/*.cpp)
file(GLOB LOCAL_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/local/*.cpp)
file(GLOB ENV_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/deconv/env/*.cpp)

# classes for topmg
file(GLOB TOPMG_FILTER_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/massmatch/*.cpp
	                     ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/oneptm/*.cpp
			     ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/diag/*.cpp)
# mass graph for topmg
file(GLOB GRAPH_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/graph/*.cpp)
file(GLOB GRAPH_ALIGN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/graphalign/*.cpp)

# MCMC for topmg
file(GLOB MCMC_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/mcmc/*.cpp)

# Deconvoluation
file(GLOB_RECURSE DECONV_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/deconv/*.cpp)
file(GLOB FEATURE_DETECT_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/feature_detect/*.cpp)

# Merge
file(GLOB MERGE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/merge/*.cpp)

# topfd_main
file(GLOB TOPFD_PROC_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_process.cpp)
file(GLOB TOPFD_CONSOLE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_argument.cpp
                               ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd.cpp)

# topfd single_main
file(GLOB TOPFD_SINGLE_PROC_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_single_process.cpp)
file(GLOB TOPFD_SINGLE_CONSOLE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_single_argument.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_single.cpp)

# toppic_main
file(GLOB TOPPIC_PROC_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_process.cpp
			     ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_argument.cpp)
file(GLOB TOPPIC_CONSOLE_SRCS	${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic.cpp)

# topmg_main
file(GLOB TOPMG_PROC_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_process.cpp
			    ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_argument.cpp)
file(GLOB TOPMG_CONSOLE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg.cpp)

# merge main
file(GLOB TOPMERGE_PROC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmerge_process.cpp
			     ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmerge_argument.cpp)
file(GLOB TOPMERGE_CONSOLE_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmerge.cpp)

# GUI
file(GLOB TOPFD_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topfd/*.cpp
                           ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topfd/*.qrc)
file(GLOB TOPPIC_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/toppic/*.cpp
			    ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/toppic/*.qrc)
file(GLOB TOPMG_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmg/*.cpp
			   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmg/*.qrc)
file(GLOB TOPMERGE_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmerge/*.cpp
			      ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmerge/*.qrc)

# unit test using catch
file(GLOB UNIT_TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/*.cpp)

add_executable(topfd 
	${TOPFD_CONSOLE_SRCS} ${TOPFD_PROC_SRCS} ${DECONV_SRCS} ${SPEC_SRCS} 
  ${FEATURE_SRCS} ${PWIZ_SRCS} ${COMMON_SRCS} ${FEATURE_DETECT_SRCS})

add_executable(topfd_single 
	${TOPFD_SINGLE_CONSOLE_SRCS} ${TOPFD_SINGLE_PROC_SRCS} ${DECONV_SRCS} ${SPEC_SRCS} 
	${PWIZ_SRCS} ${COMMON_SRCS}) 

add_executable(topfd_gui 
	${TOPFD_GUI_SRCS} ${TOPFD_PROC_SRCS} ${DECONV_SRCS} ${SPEC_SRCS} 
  ${FEATURE_SRCS} ${PWIZ_SRCS} ${COMMON_SRCS} ${FEATURE_DETECT_SRCS})

add_executable(toppic 
	${TOPPIC_CONSOLE_SRCS} ${TOPPIC_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${FILTER_SRCS}  
	${ZERO_PTM_SRCS} ${ONE_PTM_SEARCH_SRCS} ${PTM_SRCS} 
  ${TDGF_SRCS} ${LOCAL_SRCS} ${PRSM_VIEW_SRCS} ${ENV_SRCS})

add_executable(toppic_gui 
	${TOPPIC_GUI_SRCS} ${TOPPIC_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${FILTER_SRCS} ${ZERO_PTM_SRCS} 
	${ONE_PTM_SEARCH_SRCS} ${PTM_SRCS} 
  ${TDGF_SRCS} ${LOCAL_SRCS} ${PRSM_VIEW_SRCS} ${ENV_SRCS})

add_executable(toppic_test 
  ${TEST_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${FILTER_SRCS}  
	${ZERO_PTM_SRCS} ${ONE_PTM_SEARCH_SRCS} ${PTM_SRCS} 
  ${TDGF_SRCS} ${LOCAL_SRCS} ${PRSM_VIEW_SRCS} ${ENV_SRCS})

add_executable(topmg 
	${TOPMG_CONSOLE_SRCS} ${TOPMG_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${TOPMG_FILTER_SRCS} 
	${ONE_PTM_SEARCH_SRCS} ${GRAPH_SRCS} ${GRAPH_ALIGN_SRCS} ${TDGF_SRCS} 
  ${MCMC_SRCS} ${PRSM_VIEW_SRCS} ${ENV_SRCS})

add_executable(topmg_gui 
	${TOPMG_GUI_SRCS} ${TOPMG_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS}
	${TOPMG_FILTER_SRCS} ${ONE_PTM_SEARCH_SRCS} 
	${GRAPH_SRCS} ${GRAPH_ALIGN_SRCS} ${TDGF_SRCS} 
	${MCMC_SRCS} ${PRSM_VIEW_SRCS})

add_executable(topmerge 
	${TOPMERGE_CONSOLE_SRCS} ${TOPMERGE_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS}
  ${SEQ_SRCS} ${SPEC_SRCS} ${PRSM_SRCS} ${FEATURE_SRCS} ${MERGE_SRCS}
  ${ENV_SRCS})

#add_executable(topmerge_gui ${TOPMERGE_GUI_SRCS} ${TOPMERGE_GUI_QRC} ${TOPMERGE_PROCESS_SRCS} ${TOPMERGE_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS})

#add_executable(toppic_match_peak ${TOPPIC_MATCH_PEAK_SRCS} ${TOPPIC_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS} ${PRSM_VIEW_SRCS})

#add_executable(topfd_train ${TOPFD_TRAIN_SRCS} ${TOPFD_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS} ${PWIZ_SRC})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries(topfd xerces-c
        boost_program_options-mt boost_filesystem-mt
        boost_system-mt boost_thread-mt
        boost_iostreams-mt boost_chrono-mt pthread z Ws2_32)

    target_link_libraries(topfd_single xerces-c
        boost_program_options-mt boost_filesystem-mt
        boost_system-mt boost_thread-mt
        boost_iostreams-mt boost_chrono-mt pthread z Ws2_32)

    target_link_libraries(toppic xalan-c xalanMsg xerces-c
        boost_filesystem-mt boost_system-mt 
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

    target_link_libraries(topmg xalan-c xalanMsg xerces-c
        boost_filesystem-mt boost_system-mt
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

   target_link_libraries(topmerge xalan-c xerces-c boost_program_options-mt
        boost_filesystem-mt boost_system-mt boost_thread-mt pthread z Ws2_32)

    # for Qt
    SET(QTOPTION -Wl,-subsystem,windows)
    target_link_libraries(topfd_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui
        xerces-c boost_filesystem-mt boost_system-mt boost_thread-mt 
		boost_iostreams-mt boost_chrono-mt pthread z Ws2_32)

    target_link_libraries(toppic_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui xalan-c xalanMsg xerces-c
        boost_filesystem-mt boost_system-mt 
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

    target_link_libraries(topmg_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui xalan-c xalanMsg xerces-c
        boost_filesystem-mt boost_system-mt
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

#    target_link_libraries(topmerge_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui xalan-c xalanMsg xerces-c
#        boost_filesystem-mt boost_system-mt 
#        boost_program_options-mt boost_thread-mt pthread z Ws2_32)


ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows") 

# xerces-c boost_system, boost_thread boost_iostreams boost_chrono pthread z are required by
# ProteoWizard
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_link_libraries(topfd xerces-c boost_program_options
        boost_filesystem boost_system boost_thread boost_iostreams boost_chrono
        pthread z)

    target_link_libraries(topfd_gui Qt5Widgets Qt5Core Qt5Gui
        xerces-c boost_filesystem boost_system boost_thread boost_iostreams boost_chrono pthread z)

    target_link_libraries(topfd_single xerces-c boost_program_options
        boost_filesystem boost_system boost_thread boost_iostreams boost_chrono
        pthread z)

    target_link_libraries(toppic xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(toppic_gui Qt5Widgets Qt5Core Qt5Gui xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(toppic_test xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topmg xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topmg_gui Qt5Widgets Qt5Core Qt5Gui xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topmerge xalan-c xerces-c boost_program_options
	boost_filesystem boost_system boost_thread pthread z)

#    target_link_libraries(topmerge_gui Qt5Widgets Qt5Core Qt5Gui
#        xalan-c xerces-c boost_program_options boost_filesystem boost_system boost_thread boost_iostreams pthread z)

#    target_link_libraries(toppic_match_peak xalan-c xerces-c boost_program_options
#        boost_filesystem boost_system boost_thread pthread z)

#    target_link_libraries(topfd_train xalan-c xerces-c boost_program_options
#        boost_filesystem boost_system boost_thread boost_iostreams boost_chrono pthread z)

#file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/usr/bin)
install (TARGETS topfd toppic topmg topmerge DESTINATION ${CMAKE_INSTALL_PREFIX}/usr/bin/)
install (DIRECTORY toppic_resources/ DESTINATION ${CMAKE_INSTALL_PREFIX}/etc/toppic)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

