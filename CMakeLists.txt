#specify the version being used as well as the language
cmake_minimum_required(VERSION 3.5)
#Name your project here
project(TopPIC)

SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

message(STATUS "System name: " ${CMAKE_SYSTEM_NAME})

# -DWITHOUT_MZ5 is used in pwiz
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -std=c++17 -O3 -DWITHOUT_MZ5 -fPIC")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Linux specific code 
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-rpath=$ORIGIN/../lib/toppic,--disable-new-dtags -DTOPPIC_LINUX")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # Windows specific code
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTOPPIC_WINDOWS -DWIN32")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

message(STATUS "Cmake cxx flags: " ${CMAKE_CXX_FLAGS})

# for Qt
set(CMAKE_AUTOMOC ON)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
find_package(Qt5Widgets)
find_package(Qt5Core)
find_package(Qt5Gui)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext
                    ${Qt5Widgets_INCLUDE_DIRS}
                    ${Qt5Core_INCLUDE_DIRS}
                    ${Qt5Gui_INCLUDE_DIRS})

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/toppic)

# pwiz
file(GLOB PWIZ_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/utility/minimxml/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/utility/misc/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/data/common/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/data/msdata/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/pwiz/data/msdata/mzmlb/*.cpp)

# hts
file(GLOB HTS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib/*.c)

# nowide used for the Windows version of topfd
file(GLOB NOWIDE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/ext/boost/nowide/*.c)

# common
file(GLOB_RECURSE COMMON_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.cpp)
file(GLOB SEQ_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/seq/*.cpp)
file(GLOB PARA_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/para/*.cpp)
file(GLOB SPEC_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/spec/*.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/msmap/*.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/factory/*.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/env/*.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/feature/*.cpp)
file(GLOB MZML_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/ms/mzml/*.cpp)
file(GLOB PRSM_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/prsm/*.cpp)

#gui to run exe in command line
file(GLOB GUI_UTIL_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/util/*.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/common/util/*.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/common/base/*.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/common/xml/*.cpp)

# topfd 
file(GLOB_RECURSE TOPFD_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/topfd/*.cpp)

# topfd console 
file(GLOB TOPFD_CONSOLE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd_argument.cpp
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topfd.cpp)

# topfd gui
file(GLOB TOPFD_GUI_SRCS    ${CMAKE_CURRENT_SOURCE_DIR}/src/topfd/common/topfd_para.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topfd/*.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topfd/*.qrc)

# topindex 
file(GLOB FILTER_INDEX_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/massmatch/*.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/mng/*.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/index/*.cpp)

# topindex console
file(GLOB TOPINDEX_PROC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topindex_process.cpp
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topindex_argument.cpp)					 

file(GLOB TOPINDEX_CONSOLE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topindex.cpp)

# topindex gui  
file(GLOB TOPINDEX_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topindex/*.cpp
			                        ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topindex/*.qrc
                              ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topindex_argument.cpp)

# toppic
file(GLOB_RECURSE FILTER_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/filter/*.cpp)

file(GLOB SEARCH_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/search/diag/*.cpp 
                       ${CMAKE_CURRENT_SOURCE_DIR}/src/search/zeroptmsearch/*.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/src/search/oneptmsearch/*.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/src/search/varptmsearch/*.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/src/search/ptmsearch/*.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/src/search/duplicatematch/*.cpp)

file(GLOB TDGF_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/stat/count/*.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/stat/tdgf/*.cpp)
file(GLOB LOCAL_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/stat/local/*.cpp)
file(GLOB VISUAL_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/visual/*.cpp)

# toppic console
file(GLOB TOPPIC_PROC_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_process.cpp
		                  	    ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_argument.cpp
		                  	    ${CMAKE_CURRENT_SOURCE_DIR}/src/console/console_util.cpp)

file(GLOB TOPPIC_CONSOLE_SRCS	${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic.cpp)

# toppic gui 
file(GLOB TOPPIC_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/toppic/*.cpp
	                  		    ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/toppic/*.qrc
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/console/toppic_argument.cpp
		                  	    ${CMAKE_CURRENT_SOURCE_DIR}/src/console/console_util.cpp)

# topmg
file(GLOB GRAPH_ALIGN_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/search/diag/*.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/search/graph/*.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/search/graphalign/*.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/search/duplicatematch/*.cpp)

file(GLOB MCMC_SRCS  ${CMAKE_CURRENT_SOURCE_DIR}/src/stat/count/*.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/stat/mcmc/*.cpp)

# topmg console 
file(GLOB TOPMG_PROC_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_process.cpp
	                  		    ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_argument.cpp
                      			${CMAKE_CURRENT_SOURCE_DIR}/src/console/topindex_file_name.cpp
		                  	    ${CMAKE_CURRENT_SOURCE_DIR}/src/console/console_util.cpp)

file(GLOB TOPMG_CONSOLE_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg.cpp)

# topmg gui 
file(GLOB TOPMG_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmg/*.cpp
	                  		   ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topmg/*.qrc
                           ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topmg_argument.cpp
		                  	   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/console_util.cpp)

# topdiff 
file(GLOB DIFF_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/merge/*.cpp)

file(GLOB DIFF_SEARCH_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/search/duplicatematch/*.cpp)

# topdiff console 
file(GLOB TOPDIFF_PROC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff_process.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff_argument.cpp)
file(GLOB TOPDIFF_CONSOLE_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff.cpp)

# topdiff gui  
file(GLOB TOPDIFF_GUI_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/src/console/topdiff_argument.cpp
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topdiff/*.cpp
		                         ${CMAKE_CURRENT_SOURCE_DIR}/src/gui/topdiff/*.qrc)

# unit test using catch
file(GLOB UNIT_TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/test/*.cpp)


add_executable(topfd 
  ${TOPFD_CONSOLE_SRCS} ${TOPFD_SRCS} ${COMMON_SRCS} ${PARA_SRCS} 
  ${SPEC_SRCS} ${MZML_SRCS} ${PWIZ_SRCS} ${NOWIDE_SRCS})

add_executable(topindex 
  ${TOPINDEX_CONSOLE_SRCS} ${TOPINDEX_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS}
  ${SEQ_SRCS} ${PARA_SRCS} ${FILTER_INDEX_SRCS}) 

add_executable(toppic 
	${TOPPIC_CONSOLE_SRCS} ${TOPPIC_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${SEQ_SRCS} ${PARA_SRCS} ${SPEC_SRCS}  
  ${PRSM_SRCS} ${FILTER_SRCS} ${SEARCH_SRCS}  
  ${TDGF_SRCS} ${LOCAL_SRCS} ${VISUAL_SRCS})

add_executable(topmg 
	${TOPMG_CONSOLE_SRCS} ${TOPMG_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS}
  ${SEQ_SRCS} ${PARA_SRCS} ${SPEC_SRCS} ${PRSM_SRCS} ${FILTER_SRCS} 
	${GRAPH_ALIGN_SRCS} ${MCMC_SRCS} ${VISUAL_SRCS})

add_executable(topdiff 
  ${TOPDIFF_CONSOLE_SRCS} ${TOPDIFF_PROC_SRCS} ${HTS_SRCS} ${COMMON_SRCS} 
  ${DIFF_SEARCH_SRCS} ${SEQ_SRCS} ${PARA_SRCS} ${SPEC_SRCS} ${PRSM_SRCS} 
  ${DIFF_SRCS})

add_executable(topfd_gui ${TOPFD_GUI_SRCS} ${GUI_UTIL_SRCS})

add_executable(topindex_gui ${TOPINDEX_GUI_SRCS} ${GUI_UTIL_SRCS})

add_executable(toppic_gui ${TOPPIC_GUI_SRCS} ${GUI_UTIL_SRCS})

add_executable(topmg_gui ${TOPMG_GUI_SRCS} ${GUI_UTIL_SRCS})

add_executable(topdiff_gui ${TOPDIFF_GUI_SRCS} ${GUI_UTIL_SRCS})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

  add_library(onnx_lib SHARED IMPORTED)
  set_property(TARGET onnx_lib PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/windows/onnxruntime.dll)
  set_property(TARGET onnx_lib PROPERTY IMPORTED_IMPLIB ${CMAKE_CURRENT_SOURCE_DIR}/lib/windows/onnxruntime.lib)

  target_link_libraries(topfd boost_chrono-mt boost_filesystem-mt
    boost_program_options-mt boost_iostreams-mt boost_thread-mt
    xerces-c z psapi Ws2_32 onnx_lib)

  target_link_libraries(topindex boost_program_options-mt
    boost_serialization-mt boost_thread-mt xerces-c z Ws2_32)

  target_link_libraries(toppic boost_program_options-mt
    boost_serialization-mt boost_thread-mt xerces-c z Ws2_32)

  target_link_libraries(topmg boost_program_options-mt
    boost_serialization-mt boost_thread-mt xerces-c z Ws2_32)

  target_link_libraries(topdiff boost_program_options-mt
    boost_thread-mt xerces-c z Ws2_32)

  # for Qt
  SET(QTOPTION -Wl,-subsystem,windows)
  target_link_libraries(topfd_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui
    boost_thread-mt xerces-c z Ws2_32)

  target_link_libraries(topindex_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui 
    boost_program_options-mt boost_thread-mt xerces-c z Ws2_32)

  target_link_libraries(toppic_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui 
    boost_program_options-mt boost_thread-mt xerces-c z Ws2_32)

  target_link_libraries(topmg_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui 
    boost_program_options-mt boost_thread-mt xerces-c z Ws2_32)

  target_link_libraries(topdiff_gui ${QTOPTION} Qt5Widgets Qt5Core Qt5Gui 
    boost_program_options-mt boost_thread-mt xerces-c z Ws2_32)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows") 

# boost_chrono boost_file_system boost_iostreams boost_thread z are required by ProteoWizard
# boost_serialization boost_thread xerces-c z are required by toppic suite
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  message("Linux system") 
	target_link_libraries(topfd boost_chrono boost_filesystem boost_program_options 
    boost_iostreams boost_thread xerces-c z onnxruntime)

  target_link_libraries(topindex boost_program_options 
    boost_serialization boost_thread xerces-c z)

  target_link_libraries(toppic boost_program_options
    boost_serialization boost_thread xerces-c z)

  target_link_libraries(topmg boost_program_options
    boost_serialization boost_thread xerces-c z)

  target_link_libraries(topdiff boost_program_options
    boost_thread xerces-c z)

  target_link_libraries(topfd_gui Qt5Widgets Qt5Core Qt5Gui boost_thread xerces-c)

  target_link_libraries(topindex_gui Qt5Widgets Qt5Core Qt5Gui
    boost_program_options boost_thread xerces-c)

  target_link_libraries(toppic_gui Qt5Widgets Qt5Core Qt5Gui 
    boost_program_options boost_thread xerces-c)

  target_link_libraries(topmg_gui Qt5Widgets Qt5Core Qt5Gui 
    boost_program_options boost_thread xerces-c)

  target_link_libraries(topdiff_gui Qt5Widgets Qt5Core Qt5Gui
    boost_program_options boost_thread xerces-c)

install (TARGETS topfd topindex toppic topmg topdiff topfd_gui topindex_gui
  toppic_gui topmg_gui topdiff_gui DESTINATION bin/)
install (DIRECTORY resources/ DESTINATION etc/toppic)
install (DIRECTORY lib/toppic DESTINATION lib/toppic)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux") 

