#Specify the version being used as well as the language
cmake_minimum_required(VERSION 3.1)
#Name your project here
project(TopPIC)

SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../lib)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O3 -fPIC")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O2 -pg -fPIC")
message(STATUS "system name " ${CMAKE_SYSTEM_NAME})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Linux specific code
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPROT_LINUX")

  set(CMAKE_AUTOMOC ON)
  find_package(Qt5Widgets)  
  find_package(Qt5Core)
  find_package(Qt5Gui)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
                      /usr/include/libpwiz
                      /usr/include/libpwiz/libraries/boost_aux
                      ${Qt5Widgets_INCLUDE_DIRS}
                      ${Qt5Core_INCLUDE_DIRS}
                      ${Qt5Gui_INCLUDE_DIRS})
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # Windows specific code
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/windows_include)
  link_directories(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/windows_lib/)
  SET(CMAKE_LINK_LIBRARY_SUFFIXES .a)
  SET(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows") 


message(STATUS "cmake cxx flags " ${CMAKE_CXX_FLAGS})

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
file(GLOB POISSON_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/poisson/*.cpp)
file(GLOB TDGF_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/tdgf/*.cpp)
file(GLOB PRSM_VIEW_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/prsmview/*.cpp)
file(GLOB LOCAL_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/local/*.cpp)
file(GLOB TOPPIC_ARGU_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_argument.cpp)
file(GLOB TOPPIC_MAIN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic.cpp)
file(GLOB GRAPH_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/graph/*.cpp)
file(GLOB GRAPH_ALIGN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/graphalign/*.cpp)
file(GLOB TOPMG_ARGU_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_argument.cpp)
file(GLOB TOPMG_MAIN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg.cpp)
file(GLOB FEATURE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/feature/*.cpp)
file(GLOB TOPFD_ARGU_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_argument.cpp)
file(GLOB TOPFD_MAIN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd.cpp)
file(GLOB TOPFD_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topfd/*.cpp)

file(GLOB UNIT_TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/*.cpp) 

add_executable(toppic ${TOPPIC_MAIN_SRCS} ${TOPPIC_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
    ${ZERO_PTM_SRCS} ${ZERO_PTM_FILTER} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${ONE_PTM_SEARCH_SRCS}
    ${PTM_SRCS} ${TDGF_SRCS} ${LOCAL_SRCS} ${PRSM_VIEW_SRCS})

add_executable(topmg ${TOPMG_MAIN_SRC} ${TOPMG_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
    ${ZERO_PTM_FILTER} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${ZERO_PTM_SRCS} ${GRAPH_SRCS}
    ${GRAPH_ALIGN_SRCS} ${TDGF_SRCS} ${ONE_PTM_SEARCH_SRCS} ${PRSM_VIEW_SRCS})

add_executable(topfd ${TOPFD_MAIN_SRCS} ${TOPFD_ARGU_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS})

add_executable(topfd_gui ${TOPFD_GUI_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS})

add_library(topfd-dev SHARED ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS})

#add_executable(toppic-test ${UNIT_TEST_SRCS} ${HTS_SRCS} ${BASE_SRCS})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Windows specific code
    target_link_libraries(toppic xalan-c xalanMsg xerces-c
        boost_filesystem-mgw48-mt-s-1_60 boost_system-mgw48-mt-s-1_60
        boost_program_options-mgw48-mt-s-1_60 boost_thread-mgw48-mt-s-1_60 pthread z Ws2_32)

    target_link_libraries(topmg xalan-c xalanMsg xerces-c
        boost_filesystem-mgw48-mt-s-1_60 boost_system-mgw48-mt-s-1_60
        boost_program_options-mgw48-mt-s-1_60 boost_thread-mgw48-mt-s-1_60 pthread z Ws2_32)

    target_link_libraries(topfd pwiz xalan-c xalanMsg xerces-c
		boost_program_options-mgw48-mt-s-1_60 boost_filesystem-mgw48-mt-s-1_60
		boost_system-mgw48-mt-s-1_60 boost_thread-mgw48-mt-s-1_60
		boost_iostreams-mgw48-mt-s-1_60 boost_regex-mgw48-mt-s-1_60 pthread z Ws2_32)
		
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Linux specific code

    target_link_libraries(toppic xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topmg xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topfd pwiz xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

    target_link_libraries(topfd-dev pwiz xalan-c xerces-c boost_filesystem
        boost_system boost_thread pthread z)

    target_link_libraries(topfd_gui Qt5Widgets Qt5Core Qt5Gui GL pwiz
        xalan-c xerces-c boost_filesystem boost_system boost_thread pthread z)

    #target_link_libraries(toppic-test xalan-c xerces-c boost_filesystem
    #    boost_system pthread z)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

