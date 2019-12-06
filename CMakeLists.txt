#Specify the version being used as well as the language
cmake_minimum_required(VERSION 3.1)
#Name your project here
project(TopPIC)

SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../lib)

SET(CMAKE_INSTALL_PREFIX ${HOME}/toppic)

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
                      ${CMAKE_CURRENT_SOURCE_DIR}/ext
                      ${Qt5Widgets_INCLUDE_DIRS}
                      ${Qt5Core_INCLUDE_DIRS}
                      ${Qt5Gui_INCLUDE_DIRS})
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # Windows specific code
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTOPPIC_WINDOWS -DWIN32")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
                      ${CMAKE_CURRENT_SOURCE_DIR}/ext
                      ${Qt5Widgets_INCLUDE_DIRS}
                      ${Qt5Core_INCLUDE_DIRS}
                      ${Qt5Gui_INCLUDE_DIRS})
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

message(STATUS "cmake cxx flags " ${CMAKE_CXX_FLAGS})

# pwiz
file(GLOB PWIZ_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/utility/minimxml/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/utility/misc/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/data/common/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/data/msdata/*.cpp)

# hts
file(GLOB HTS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib/*.c)

# common
file(GLOB_RECURSE COMMON_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp)
file(GLOB SEQ_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/seq/*.cpp)
file(GLOB SPEC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/spec/*.cpp)
file(GLOB ENV_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/env/*.cpp)
file(GLOB FEATURE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/feature/*.cpp)
file(GLOB PRSM_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/prsm/*.cpp)

# classes for topfd 
file(GLOB_RECURSE TOPFD_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/topfd/*.cpp)

# topfd console 
file(GLOB TOPFD_CONSOLE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_argument.cpp
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd.cpp)
# topfd gui
file(GLOB TOPFD_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topfd/*.cpp
                           ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topfd/*.qrc)

# topfd_single console main
file(GLOB TOPFD_SINGLE_PROC_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_single_process.cpp)
file(GLOB TOPFD_SINGLE_CONSOLE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_single_argument.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_single.cpp)

# classes for toppic
file(GLOB_RECURSE FILTER_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/*.cpp)

file(GLOB ZERO_PTM_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/search/zeroptmsearch/*.cpp)
file(GLOB ONE_PTM_SEARCH_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/search/oneptmsearch/*.cpp)
file(GLOB PTM_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/search/ptmsearch/*.cpp)

file(GLOB TDGF_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/stat/tdgf/*.cpp)
file(GLOB LOCAL_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/stat/local/*.cpp)
file(GLOB VISUAL_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/visual/*.cpp)

# toppic console main
file(GLOB TOPPIC_PROC_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_process.cpp
			     ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_argument.cpp)
file(GLOB TOPPIC_CONSOLE_SRCS	${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic.cpp)

# toppic gui main
file(GLOB TOPPIC_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/toppic/*.cpp
			    ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/toppic/*.qrc)

# classes for topmg
file(GLOB TOPMG_FILTER_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/massmatch/*.cpp
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/oneptm/*.cpp
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/diag/*.cpp)

file(GLOB GRAPH_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/search/graph/*.cpp)
file(GLOB GRAPH_ALIGN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/search/graphalign/*.cpp)

file(GLOB MCMC_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/stat/mcmc/*.cpp)

# topmg console main
file(GLOB TOPMG_PROC_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_process.cpp
			    ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_argument.cpp)
file(GLOB TOPMG_CONSOLE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg.cpp)

# topmg gui main
file(GLOB TOPMG_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmg/*.cpp
			   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmg/*.qrc)

# topdiff classes
file(GLOB DIFF_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/merge/*.cpp)

# topdiff console main
file(GLOB TOPDIFF_PROC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff_process.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff_argument.cpp)
file(GLOB TOPDIFF_CONSOLE_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff.cpp)

# topdiff gui main 
file(GLOB TOPDIFF_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topdiff/*.cpp
			      ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topdiff/*.qrc)

# unit test using catch
file(GLOB UNIT_TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/*.cpp)

add_executable(topfd 
  ${TOPFD_CONSOLE_SRCS} ${TOPFD_SRCS} ${SPEC_SRCS} ${ENV_SRCS} 
  ${FEATURE_SRCS} ${PWIZ_SRCS} ${COMMON_SRCS})

add_executable(topfd_single 
  ${TOPFD_SINGLE_CONSOLE_SRCS} ${TOPFD_SINGLE_PROC_SRCS} ${TOPFD_SRCS}
  ${SPEC_SRCS} ${ENV_SRCS} ${PWIZ_SRCS} ${COMMON_SRCS}) 

add_executable(topfd_gui 
  ${TOPFD_GUI_SRCS} ${TOPFD_SRCS} ${SPEC_SRCS} ${ENV_SRCS} 
  ${FEATURE_SRCS} ${PWIZ_SRCS} ${COMMON_SRCS})

add_executable(toppic 
	${TOPPIC_CONSOLE_SRCS} ${TOPPIC_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${ENV_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${FILTER_SRCS}  
	${ZERO_PTM_SRCS} ${ONE_PTM_SEARCH_SRCS} ${PTM_SRCS} 
  ${TDGF_SRCS} ${LOCAL_SRCS} ${VISUAL_SRCS})

add_executable(toppic_gui 
	${TOPPIC_GUI_SRCS} ${TOPPIC_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${ENV_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${FILTER_SRCS} 
  ${ZERO_PTM_SRCS} ${ONE_PTM_SEARCH_SRCS} ${PTM_SRCS} 
  ${TDGF_SRCS} ${LOCAL_SRCS} ${VISUAL_SRCS})


add_executable(topmg 
	${TOPMG_CONSOLE_SRCS} ${TOPMG_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${ENV_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${TOPMG_FILTER_SRCS} 
	${ONE_PTM_SEARCH_SRCS} ${GRAPH_SRCS} ${GRAPH_ALIGN_SRCS} ${TDGF_SRCS} 
  ${MCMC_SRCS} ${VISUAL_SRCS})

add_executable(topmg_gui 
	${TOPMG_GUI_SRCS} ${TOPMG_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${ENV_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS}
	${TOPMG_FILTER_SRCS} ${ONE_PTM_SEARCH_SRCS} 
	${GRAPH_SRCS} ${GRAPH_ALIGN_SRCS} ${TDGF_SRCS} 
  ${MCMC_SRCS} ${VISUAL_SRCS})

add_executable(topdiff 
  ${TOPDIFF_CONSOLE_SRCS} ${TOPDIFF_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS}
  ${SEQ_SRCS} ${SPEC_SRCS} ${ENV_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${DIFF_SRCS})

add_executable(topdiff_gui ${TOPDIFF_GUI_SRCS} ${TOPDIFF_GUI_QRC}
  ${TOPDIFF_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${SPEC_SRCS} ${ENV_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${DIFF_SRCS})

#add_executable(toppic_test 
#  ${TEST_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
#  ${SEQ_SRCS} ${SPEC_SRCS} ${ENV_SRCS} ${FEATURE_SRCS} ${PRSM_SRCS} ${FILTER_SRCS}  
#  ${ZERO_PTM_SRCS} ${ONE_PTM_SEARCH_SRCS} ${PTM_SRCS} 
#  ${TDGF_SRCS} ${LOCAL_SRCS} ${VISUAL_SRCS})


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

    target_link_libraries(toppic xerces-c
        boost_filesystem-mt boost_system-mt 
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

    target_link_libraries(topmg xerces-c
        boost_filesystem-mt boost_system-mt
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

   target_link_libraries(topdiff xerces-c boost_program_options-mt
        boost_filesystem-mt boost_system-mt boost_thread-mt pthread z Ws2_32)

    # for Qt
    SET(QTOPTION -Wl,-subsystem,windows)
    target_link_libraries(topfd_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui
        xerces-c boost_filesystem-mt boost_system-mt boost_thread-mt 
	    	boost_iostreams-mt boost_chrono-mt pthread z Ws2_32)

    target_link_libraries(toppic_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui xerces-c
        boost_filesystem-mt boost_system-mt 
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

    target_link_libraries(topmg_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui xerces-c
        boost_filesystem-mt boost_system-mt
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)

    target_link_libraries(topdiff_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui xerces-c
        boost_filesystem-mt boost_system-mt 
        boost_program_options-mt boost_thread-mt pthread z Ws2_32)


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

    target_link_libraries(toppic xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(toppic_gui Qt5Widgets Qt5Core Qt5Gui xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topmg xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topmg_gui Qt5Widgets Qt5Core Qt5Gui xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topdiff xerces-c boost_program_options
	      boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topdiff_gui Qt5Widgets Qt5Core Qt5Gui
        xerces-c boost_program_options boost_filesystem boost_system boost_thread boost_iostreams pthread z)

#    target_link_libraries(toppic_test xerces-c boost_program_options
#  boost_filesystem boost_system boost_thread pthread z)

#    target_link_libraries(toppic_match_peak xerces-c boost_program_options
#        boost_filesystem boost_system boost_thread pthread z)

#    target_link_libraries(topfd_train xerces-c boost_program_options
#        boost_filesystem boost_system boost_thread boost_iostreams boost_chrono pthread z)

#file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/usr/bin)
install (TARGETS topfd toppic topmg topdiff DESTINATION ${CMAKE_INSTALL_PREFIX}/usr/bin/)
install (DIRECTORY toppic_resources/ DESTINATION ${CMAKE_INSTALL_PREFIX}/etc/toppic)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

