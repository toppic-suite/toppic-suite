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

file(GLOB BASE_SCRS ${CMAKE_CURRENT_SOURCE_DIR}/src/base/*.cpp)
file(GLOB SPEC_SCRS ${CMAKE_CURRENT_SOURCE_DIR}/src/spec/*.cpp)
file(GLOB PRSM_SCRS  ${CMAKE_CURRENT_SOURCE_DIR}/src/prsm/*.cpp)
file(GLOB ZERO_PTM_SCRS  ${CMAKE_CURRENT_SOURCE_DIR}/src/zeroptmsearch/*.cpp)
file(GLOB FILTER_SCRS  ${CMAKE_CURRENT_SOURCE_DIR}/src/filterdiagonal/*.cpp)
file(GLOB PTM_SCRS   ${CMAKE_CURRENT_SOURCE_DIR}/src/ptmsearch/*.cpp)
file(GLOB POISSON_SCRS   ${CMAKE_CURRENT_SOURCE_DIR}/src/poisson/*.cpp)
file(GLOB TDGF_SCRS   ${CMAKE_CURRENT_SOURCE_DIR}/src/tdgf/*.cpp)
file(GLOB XPP_SCRS   ${CMAKE_CURRENT_SOURCE_DIR}/src/xpp/*.cpp)
file(GLOB TEST_SCRS   ${CMAKE_CURRENT_SOURCE_DIR}/src/test/argument.cpp)
file(GLOB MAIN_SCRS   ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_console.cpp) 
file(GLOB MULTIPLE_SCRS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_multiple.cpp) 
file(GLOB TEST_EVALUE_SCRS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/test_evalue.cpp) 

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/linux_lib/)

message(STATUS "multiple sources " ${MULTIPLE_SCRS})

add_executable(topid ${BASE_SCRS} ${SPEC_SCRS} ${PRSM_SCRS}
        ${ZERO_PTM_SCRS} ${FILTER_SCRS} ${PTM_SCRS} ${POISSON_SCRS}
        ${TDGF_SCRS} ${XPP_SCRS} ${TEST_SCRS} ${MAIN_SCRS})

add_executable(test_multiple ${BASE_SCRS} ${SPEC_SCRS} ${PRSM_SCRS}
        ${ZERO_PTM_SCRS} ${FILTER_SCRS} ${PTM_SCRS} ${POISSON_SCRS} 
        ${TDGF_SCRS} ${XPP_SCRS} ${TEST_SCRS} ${MULTIPLE_SCRS})

add_executable(test_evalue ${BASE_SCRS} ${SPEC_SCRS} ${PRSM_SCRS}
        ${ZERO_PTM_SCRS} ${FILTER_SCRS} ${PTM_SCRS} ${POISSON_SCRS} 
        ${TDGF_SCRS} ${XPP_SCRS} ${TEST_SCRS} ${TEST_EVALUE_SCRS})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Linux specific code
  target_link_libraries(topid xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)
  target_link_libraries(test_multiple xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)
  target_link_libraries(test_evalue xalan-c xalanMsg xerces-c boost_program_options boost_filesystem boost_system pthread)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Mingw")
  # Mingw specific code
  target_link_libraries(topid xalan-c xalanMsg xerces-c
          boost_filesystem-mgw48-1_55 boost_system-mgw48-1_55
          boost_program_options-mgw48-1_55)
  target_link_libraries(test_multiple xalan-c xalanMsg xerces-c
          boost_filesystem-mgw48-1_55 boost_system-mgw48-1_55
          boost_program_options-mgw48-1_55)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Mingw") 

