#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)
#Name your project here
project(TopPIC)

SET( CMAKE_VERBOSE_MAKEFILE ON )

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O3")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O2 -pg ")
message(STATUS "system name " ${CMAKE_SYSTEM_NAME})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Linux specific code
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPROT_LINUX")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
                      /usr/include/libpwiz
                      /usr/include/libpwiz/libraries/boost_aux)
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
file(GLOB TOPFD_MAIN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/feature/topfd_cmd.cpp)
file(GLOB GRAPH_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/graph/*.cpp)
file(GLOB GRAPH_ALIGN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/graphalign/*.cpp)
file(GLOB TOPMG_ARGU_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_argument.cpp)
file(GLOB GRAPH_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg.cpp)
file(GLOB FEATURE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/feature/*.cpp)
list(REMOVE_ITEM FEATURE_SRCS ${TOPFD_MAIN_SRCS})
file(GLOB TEST_ZERO_PTM_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_zero_ptm.cpp) 
file(GLOB TEST_DIAG_FILTER_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_diag_filter.cpp) 
file(GLOB TEST_PTM_SEARCH_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_ptm_search.cpp) 
file(GLOB TEST_EVALUE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_evalue.cpp) 
file(GLOB TEST_PRSM_VIEW_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_prsm_view.cpp) 
file(GLOB TEST_BSE_OPT ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_base_opt.cpp) 
file(GLOB TEST_GRAPH ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_graph_align.cpp) 

add_executable(toppic ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
    ${ZERO_PTM_SRCS} ${ZERO_PTM_FILTER} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${ONE_PTM_SEARCH_SRCS}
    ${PTM_SRCS} ${TDGF_SRCS} ${LOCAL_SRCS} ${PRSM_VIEW_SRCS} ${TOPPIC_ARGU_SRCS} ${TOPPIC_MAIN_SRCS})

add_executable(topmg ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS} ${ZERO_PTM_FILTER} ${DIAG_FILTER_SRCS}
    ${ONE_PTM_FILTER_SRCS} ${ZERO_PTM_SRCS} ${GRAPH_SRCS} ${GRAPH_ALIGN_SRCS} ${TDGF_SRCS}
    ${TOPMG_ARGU_SRCS} ${ONE_PTM_SEARCH_SRCS} ${PRSM_VIEW_SRCS} ${GRAPH_MAIN})

add_executable(topfd ${TOPFD_MAIN_SRCS} ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS})

add_library(topfd-dev SHARED ${HTS_SRCS} ${BASE_SRCS} ${SPEC_SRCS} ${FEATURE_SRCS})

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

    target_link_libraries(topfd-dev pwiz xalan-c xerces-c boost_program_options
        boost_filesystem boost_system boost_thread pthread z)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

