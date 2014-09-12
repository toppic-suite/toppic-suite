#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)
#Name your project here
project(TopID)

SET( CMAKE_VERBOSE_MAKEFILE ON )

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)


SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O3 ")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pg ")
message(STATUS "system name " ${CMAKE_SYSTEM_NAME})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Linux specific code
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPROT_LINUX")
  SET(CMAKE_EXE_LINKER_FLAGS  "-static")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/linux_include)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Mingw")
  # Mingw specific code
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/windows_include)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Mingw") 


message(STATUS "cmake cxx flags " ${CMAKE_CXX_FLAGS})

file(GLOB BASE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/base/*.cpp)
file(GLOB SPEC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/spec/*.cpp)
file(GLOB PRSM_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/prsm/*.cpp)
file(GLOB ZERO_PTM_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/zeroptmsearch/*.cpp)
file(GLOB DIAG_FILTER_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/diagfilter/*.cpp)
file(GLOB ONE_PTM_FILTER_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/oneptmfilter/*.cpp)
file(GLOB PTM_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/ptmsearch/*.cpp)
file(GLOB POISSON_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/poisson/*.cpp)
file(GLOB TDGF_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/tdgf/*.cpp)
file(GLOB PRSM_VIEW_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/prsmview/*.cpp)
file(GLOB ARGU_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/argument.cpp)
file(GLOB MAIN_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic.cpp) 
file(GLOB TEST_ZERO_PTM_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_zero_ptm.cpp) 
file(GLOB TEST_DIAG_FILTER_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_diag_filter.cpp) 
file(GLOB TEST_PTM_SEARCH_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_ptm_search.cpp) 
file(GLOB TEST_EVALUE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_evalue.cpp) 
file(GLOB TEST_PRSM_VIEW_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_prsm_view.cpp) 

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/linux_lib/)

add_executable(toppic ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
        ${ZERO_PTM_SRCS} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${PTM_SRCS} ${POISSON_SRCS}
        ${TDGF_SRCS} ${PRSM_VIEW_SRCS} ${ARGU_SRCS} ${MAIN_SRCS})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Mingw")
  # Mingw specific code
  target_link_libraries(toppc xalan-c xalanMsg xerces-c
          boost_filesystem-mgw48-1_55 boost_system-mgw48-1_55
          boost_program_options-mgw48-1_55)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Mingw") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Linux specific code
  target_link_libraries(toppic xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)

  add_executable(test_zero_ptm ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
          ${ZERO_PTM_SRCS} ${ARGU_SRCS} ${TEST_ZERO_PTM_SRCS})
  target_link_libraries(test_zero_ptm xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)

  add_executable(test_diag_filter ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
          ${DIAG_FILTER_SRCS} ${ARGU_SRCS} ${TEST_DIAG_FILTER_SRCS})
  target_link_libraries(test_diag_filter xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)

  add_executable(test_ptm_search ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
          ${ZERO_PTM_SRCS} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS}
          ${PTM_SRCS} ${ARGU_SRCS} ${TEST_PTM_SEARCH_SRCS})
  target_link_libraries(test_ptm_search xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)

  add_executable(test_evalue ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
          ${ZERO_PTM_SRCS} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${PTM_SRCS}  
          ${TDGF_SRCS} ${ARGU_SRCS} ${TEST_EVALUE_SRCS})
  target_link_libraries(test_evalue xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)

  add_executable(test_prsm_view ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
          ${ZERO_PTM_SRCS} ${DIAG_FILTER_SRCS} ${ONE_PTM_FILTER_SRCS} ${PTM_SRCS}  
          ${TDGF_SRCS} ${PRSM_VIEW_SRCS} ${ARGU_SRCS} ${TEST_PRSM_VIEW_SRCS})
  target_link_libraries(test_prsm_view xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)

  #  add_executable(test_multiple ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
  #        ${ZERO_PTM_SRCS} ${FILTER_SRCS} ${PTM_SRCS} ${POISSON_SRCS} 
  #        ${TDGF_SRCS} ${XPP_SRCS} ${ARGU_SRCS} ${MULTIPLE_SRCS})
  #  add_executable(test_evalue ${BASE_SRCS} ${SPEC_SRCS} ${PRSM_SRCS}
  #        ${ZERO_PTM_SRCS} ${FILTER_SRCS} ${PTM_SRCS} ${POISSON_SRCS} 
  #        ${TDGF_SRCS} ${XPP_SRCS} ${ARGU_SRCS} ${TEST_EVALUE_SRCS})
  #  target_link_libraries(test_multiple xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)
  #  target_link_libraries(test_evalue xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 


